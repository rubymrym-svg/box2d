// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT

#include "cluster.h"

#include "arena_allocator.h"
#include "array.h"
#include "body.h"
#include "constraint_graph.h"
#include "contact.h"
#include "contact_solver.h"
#include "joint.h"
#include "physics_world.h"
#include "solver.h"
#include "solver_set.h"

#include <string.h>

void b2CreateClusters( b2ClusterManager* manager )
{
	*manager = (b2ClusterManager){ 0 };

	for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
	{
		b2Cluster* cluster = manager->clusters + i;
		cluster->bodyIndices = b2IntArray_Create( 16 );
	}
}

void b2DestroyClusters( b2ClusterManager* manager )
{
	for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
	{
		b2Cluster* cluster = manager->clusters + i;
		b2IntArray_Destroy( &cluster->bodyIndices );
	}
}

void b2ComputeClusters( b2World* world )
{
	b2ClusterManager* manager = &world->clusterManager;
	b2Cluster* clusters = manager->clusters;

	b2SolverSet* awakeSet = b2SolverSetArray_Get( &world->solverSets, b2_awakeSet );
	int awakeCount = awakeSet->bodySims.count;
	b2BodySim* bodySims = awakeSet->bodySims.data;

	for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
	{
		b2IntArray_Clear( &clusters[i].bodyIndices );
	}

	if ( awakeCount == 0 )
	{
		return;
	}

	// First-time or too few bodies: full k-means with seeding
	if ( manager->initialized == false || awakeCount < B2_CLUSTER_COUNT )
	{
		int seedCount = b2MinInt( awakeCount, B2_CLUSTER_COUNT );
		for ( int i = 0; i < seedCount; ++i )
		{
			clusters[i].center = bodySims[i].center;
			b2IntArray_Push( &clusters[i].bodyIndices, i );
			bodySims[i].clusterIndex = i;
		}

		if ( awakeCount < B2_CLUSTER_COUNT )
		{
			return;
		}

		for ( int iteration = 0; iteration < 32; ++iteration )
		{
			for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
			{
				clusters[i].accumulator = b2Vec2_zero;
				clusters[i].bodyIndices.count = 0;
			}

			for ( int i = 0; i < awakeCount; ++i )
			{
				b2Vec2 p = bodySims[i].center;

				float minDistanceSquared = b2DistanceSquared( p, clusters[0].center );
				int bestIndex = 0;

				for ( int j = 1; j < B2_CLUSTER_COUNT; ++j )
				{
					float distanceSquared = b2DistanceSquared( p, clusters[j].center );
					if ( distanceSquared < minDistanceSquared )
					{
						minDistanceSquared = distanceSquared;
						bestIndex = j;
					}
				}

				bodySims[i].clusterIndex = bestIndex;
				b2IntArray_Push( &clusters[bestIndex].bodyIndices, i );
				clusters[bestIndex].accumulator = b2Add( clusters[bestIndex].accumulator, p );
			}

			for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
			{
				int clusterBodyCount = clusters[i].bodyIndices.count;
				if ( clusterBodyCount > 0 )
				{
					clusters[i].center = b2MulSV( 1.0f / clusterBodyCount, clusters[i].accumulator );
				}
			}
		}

		manager->initialized = true;
		return;
	}

	// Incremental: refine clusters using persistent centers
	for ( int iteration = 0; iteration < 4; ++iteration )
	{
		for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
		{
			clusters[i].accumulator = b2Vec2_zero;
			clusters[i].bodyIndices.count = 0;
		}

		// Assign each body to nearest cluster center
		for ( int i = 0; i < awakeCount; ++i )
		{
			b2Vec2 p = bodySims[i].center;

			float minDistanceSquared = b2DistanceSquared( p, clusters[0].center );
			int bestIndex = 0;

			for ( int j = 1; j < B2_CLUSTER_COUNT; ++j )
			{
				float distanceSquared = b2DistanceSquared( p, clusters[j].center );
				if ( distanceSquared < minDistanceSquared )
				{
					minDistanceSquared = distanceSquared;
					bestIndex = j;
				}
			}

			bodySims[i].clusterIndex = bestIndex;
			b2IntArray_Push( &clusters[bestIndex].bodyIndices, i );
			clusters[bestIndex].accumulator = b2Add( clusters[bestIndex].accumulator, p );
		}

		// Update centers, handle empty clusters
		for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
		{
			int clusterBodyCount = clusters[i].bodyIndices.count;
			if ( clusterBodyCount > 0 )
			{
				clusters[i].center = b2MulSV( 1.0f / clusterBodyCount, clusters[i].accumulator );
			}
			else
			{
				// Re-seed empty cluster from body furthest from its assigned center
				float maxDistanceSquared = -1.0f;
				int maxBody = 0;
				for ( int b = 0; b < awakeCount; ++b )
				{
					int ci = bodySims[b].clusterIndex;
					float d = b2DistanceSquared( bodySims[b].center, clusters[ci].center );
					if ( d > maxDistanceSquared )
					{
						maxDistanceSquared = d;
						maxBody = b;
					}
				}
				clusters[i].center = bodySims[maxBody].center;
			}
		}
	}
}

// Convert a pair (a, b) with a < b into a linear border index
static inline int b2GetBorderIndex( int a, int b )
{
	B2_ASSERT( a < b );
	// Row-major upper triangular index: sum of (B2_CLUSTER_COUNT - 1) + ... + (B2_CLUSTER_COUNT - a) + (b - a - 1)
	return a * ( 2 * B2_CLUSTER_COUNT - a - 1 ) / 2 + ( b - a - 1 );
}

void b2ClassifyConstraints( b2World* world, b2StepContext* context )
{
	b2ConstraintGraph* graph = &world->constraintGraph;
	b2SolverSet* awakeSet = b2SolverSetArray_Get( &world->solverSets, b2_awakeSet );
	b2BodySim* bodySims = awakeSet->bodySims.data;

	// Temporary counts for sizing
	int clusterContactCounts[B2_CLUSTER_COUNT] = { 0 };
	int clusterJointCounts[B2_CLUSTER_COUNT] = { 0 };

	// Use a flat array for border counts, indexed by b2GetBorderIndex
	int borderContactCounts[B2_MAX_BORDERS] = { 0 };
	int borderJointCounts[B2_MAX_BORDERS] = { 0 };

	// First pass: count contacts per cluster/border across all graph colors
	for ( int colorIndex = 0; colorIndex < B2_GRAPH_COLOR_COUNT - 1; ++colorIndex )
	{
		b2GraphColor* color = graph->colors + colorIndex;
		int contactCount = color->contactSims.count;
		b2ContactSim* contacts = color->contactSims.data;

		for ( int i = 0; i < contactCount; ++i )
		{
			b2ContactSim* contactSim = contacts + i;
			int indexA = contactSim->bodySimIndexA;
			int indexB = contactSim->bodySimIndexB;

			// Static bodies have B2_NULL_INDEX for bodySimIndex
			if ( indexA == B2_NULL_INDEX || indexB == B2_NULL_INDEX )
			{
				// Static-dynamic contact: assign to dynamic body's cluster
				int dynamicIndex = ( indexA != B2_NULL_INDEX ) ? indexA : indexB;
				int clusterIdx = bodySims[dynamicIndex].clusterIndex;
				clusterContactCounts[clusterIdx] += 1;
			}
			else
			{
				int clusterA = bodySims[indexA].clusterIndex;
				int clusterB = bodySims[indexB].clusterIndex;

				if ( clusterA == clusterB )
				{
					clusterContactCounts[clusterA] += 1;
				}
				else
				{
					int a = clusterA < clusterB ? clusterA : clusterB;
					int b = clusterA < clusterB ? clusterB : clusterA;
					int borderIdx = b2GetBorderIndex( a, b );
					borderContactCounts[borderIdx] += 1;
				}
			}
		}
	}

	// First pass: count joints per cluster/border across all graph colors
	for ( int colorIndex = 0; colorIndex < B2_GRAPH_COLOR_COUNT - 1; ++colorIndex )
	{
		b2GraphColor* color = graph->colors + colorIndex;
		int jointCount = color->jointSims.count;
		b2JointSim* joints = color->jointSims.data;

		for ( int i = 0; i < jointCount; ++i )
		{
			b2JointSim* jointSim = joints + i;
			int bodyIdA = jointSim->bodyIdA;
			int bodyIdB = jointSim->bodyIdB;

			b2Body* bodyA = b2BodyArray_Get( &world->bodies, bodyIdA );
			b2Body* bodyB = b2BodyArray_Get( &world->bodies, bodyIdB );

			int indexA = ( bodyA->setIndex == b2_awakeSet ) ? bodyA->localIndex : B2_NULL_INDEX;
			int indexB = ( bodyB->setIndex == b2_awakeSet ) ? bodyB->localIndex : B2_NULL_INDEX;

			if ( indexA == B2_NULL_INDEX || indexB == B2_NULL_INDEX )
			{
				int dynamicIndex = ( indexA != B2_NULL_INDEX ) ? indexA : indexB;
				int clusterIdx = bodySims[dynamicIndex].clusterIndex;
				clusterJointCounts[clusterIdx] += 1;
			}
			else
			{
				int clusterA = bodySims[indexA].clusterIndex;
				int clusterB = bodySims[indexB].clusterIndex;

				if ( clusterA == clusterB )
				{
					clusterJointCounts[clusterA] += 1;
				}
				else
				{
					int a = clusterA < clusterB ? clusterA : clusterB;
					int b = clusterA < clusterB ? clusterB : clusterA;
					int borderIdx = b2GetBorderIndex( a, b );
					borderJointCounts[borderIdx] += 1;
				}
			}
		}
	}

	// Allocate cluster solve data arrays from arena
	b2ClusterManager* manager = &world->clusterManager;
	for ( int i = 0; i < B2_CLUSTER_COUNT; ++i )
	{
		b2ClusterSolveData* cd = context->clusterData + i;
		b2Cluster* cluster = manager->clusters + i;

		int cc = clusterContactCounts[i];
		int jc = clusterJointCounts[i];

		cd->contacts = ( cc > 0 )
			? b2AllocateArenaItem( &world->arena, cc * sizeof( b2ContactSim* ), "cluster contacts" )
			: NULL;
		cd->contactCount = 0;

		cd->joints = ( jc > 0 )
			? b2AllocateArenaItem( &world->arena, jc * sizeof( b2JointSim* ), "cluster joints" )
			: NULL;
		cd->jointCount = 0;

		cd->contactConstraints = ( cc > 0 )
			? b2AllocateArenaItem( &world->arena, cc * sizeof( b2ContactConstraint ), "cluster contact constraints" )
			: NULL;

		// Allocate per-cluster local body state array for L1 cache locality
		int bodyCount = cluster->bodyIndices.count;
		cd->bodyCount = bodyCount;
		cd->bodyIndices = cluster->bodyIndices.data;
		cd->localStates = ( bodyCount > 0 )
			? b2AllocateArenaItem( &world->arena, bodyCount * sizeof( b2BodyState ), "cluster local states" )
			: NULL;

		// Build reverse mapping: global awake index -> local cluster index
		for ( int k = 0; k < bodyCount; ++k )
		{
			int globalIndex = cluster->bodyIndices.data[k];
			bodySims[globalIndex].localClusterIndex = k;
		}

		b2AtomicStoreInt( &cd->solveComplete, 0 );
	}

	// Count non-empty borders and allocate
	int borderCount = 0;
	for ( int i = 0; i < B2_MAX_BORDERS; ++i )
	{
		if ( borderContactCounts[i] > 0 || borderJointCounts[i] > 0 )
		{
			borderCount += 1;
		}
	}

	context->borderCount = borderCount;
	context->borders = ( borderCount > 0 )
		? b2AllocateArenaItem( &world->arena, borderCount * sizeof( b2BorderConstraints ), "border constraints" )
		: NULL;

	// Fill in border data
	int borderWriteIndex = 0;
	for ( int a = 0; a < B2_CLUSTER_COUNT; ++a )
	{
		for ( int b = a + 1; b < B2_CLUSTER_COUNT; ++b )
		{
			int flatIdx = b2GetBorderIndex( a, b );
			int cc = borderContactCounts[flatIdx];
			int jc = borderJointCounts[flatIdx];

			if ( cc == 0 && jc == 0 )
			{
				continue;
			}

			b2BorderConstraints* border = context->borders + borderWriteIndex;
			border->clusterA = a;
			border->clusterB = b;

			border->contacts = ( cc > 0 )
				? b2AllocateArenaItem( &world->arena, cc * sizeof( b2ContactSim* ), "border contacts" )
				: NULL;
			border->contactCount = 0;

			border->joints = ( jc > 0 )
				? b2AllocateArenaItem( &world->arena, jc * sizeof( b2JointSim* ), "border joints" )
				: NULL;
			border->jointCount = 0;

			border->contactConstraints = ( cc > 0 )
				? b2AllocateArenaItem( &world->arena, cc * sizeof( b2ContactConstraint ), "border contact constraints" )
				: NULL;

			// Store the border write index in the flat array for the second pass
			// Reuse borderContactCounts as a mapping from flat index -> border write index
			borderContactCounts[flatIdx] = borderWriteIndex;
			borderWriteIndex += 1;
		}
	}

	// Second pass: distribute contact pointers to clusters and borders
	for ( int colorIndex = 0; colorIndex < B2_GRAPH_COLOR_COUNT - 1; ++colorIndex )
	{
		b2GraphColor* color = graph->colors + colorIndex;
		int contactCount = color->contactSims.count;
		b2ContactSim* contacts = color->contactSims.data;

		for ( int i = 0; i < contactCount; ++i )
		{
			b2ContactSim* contactSim = contacts + i;
			int indexA = contactSim->bodySimIndexA;
			int indexB = contactSim->bodySimIndexB;

			if ( indexA == B2_NULL_INDEX || indexB == B2_NULL_INDEX )
			{
				int dynamicIndex = ( indexA != B2_NULL_INDEX ) ? indexA : indexB;
				int clusterIdx = bodySims[dynamicIndex].clusterIndex;
				b2ClusterSolveData* cd = context->clusterData + clusterIdx;
				cd->contacts[cd->contactCount++] = contactSim;
			}
			else
			{
				int clusterA = bodySims[indexA].clusterIndex;
				int clusterB = bodySims[indexB].clusterIndex;

				if ( clusterA == clusterB )
				{
					b2ClusterSolveData* cd = context->clusterData + clusterA;
					cd->contacts[cd->contactCount++] = contactSim;
				}
				else
				{
					int a = clusterA < clusterB ? clusterA : clusterB;
					int b = clusterA < clusterB ? clusterB : clusterA;
					int flatIdx = b2GetBorderIndex( a, b );
					// borderContactCounts was repurposed as the mapping
					int bIdx = borderContactCounts[flatIdx];
					b2BorderConstraints* border = context->borders + bIdx;
					border->contacts[border->contactCount++] = contactSim;
				}
			}
		}
	}

	// Second pass: distribute joint pointers to clusters and borders
	for ( int colorIndex = 0; colorIndex < B2_GRAPH_COLOR_COUNT - 1; ++colorIndex )
	{
		b2GraphColor* color = graph->colors + colorIndex;
		int jointCount = color->jointSims.count;
		b2JointSim* joints = color->jointSims.data;

		for ( int i = 0; i < jointCount; ++i )
		{
			b2JointSim* jointSim = joints + i;
			int bodyIdA = jointSim->bodyIdA;
			int bodyIdB = jointSim->bodyIdB;

			b2Body* bodyA = b2BodyArray_Get( &world->bodies, bodyIdA );
			b2Body* bodyB = b2BodyArray_Get( &world->bodies, bodyIdB );

			int indexA = ( bodyA->setIndex == b2_awakeSet ) ? bodyA->localIndex : B2_NULL_INDEX;
			int indexB = ( bodyB->setIndex == b2_awakeSet ) ? bodyB->localIndex : B2_NULL_INDEX;

			if ( indexA == B2_NULL_INDEX || indexB == B2_NULL_INDEX )
			{
				int dynamicIndex = ( indexA != B2_NULL_INDEX ) ? indexA : indexB;
				int clusterIdx = bodySims[dynamicIndex].clusterIndex;
				b2ClusterSolveData* cd = context->clusterData + clusterIdx;
				cd->joints[cd->jointCount++] = jointSim;
			}
			else
			{
				int clusterA = bodySims[indexA].clusterIndex;
				int clusterB = bodySims[indexB].clusterIndex;

				if ( clusterA == clusterB )
				{
					b2ClusterSolveData* cd = context->clusterData + clusterA;
					cd->joints[cd->jointCount++] = jointSim;
				}
				else
				{
					int a = clusterA < clusterB ? clusterA : clusterB;
					int b = clusterA < clusterB ? clusterB : clusterA;
					int flatIdx = b2GetBorderIndex( a, b );
					int bIdx = borderContactCounts[flatIdx];
					b2BorderConstraints* border = context->borders + bIdx;
					border->joints[border->jointCount++] = jointSim;
				}
			}
		}
	}
}
