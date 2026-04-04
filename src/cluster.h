// SPDX-FileCopyrightText: 2025 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "array.h"
#include "atomic.h"
#include "box2d/constants.h"
#include "box2d/math_functions.h"

typedef struct b2ContactConstraint b2ContactConstraint;
typedef struct b2ContactSim b2ContactSim;
typedef struct b2JointSim b2JointSim;
typedef struct b2StepContext b2StepContext;
typedef struct b2World b2World;

// Maximum number of cluster pair borders: C(16,2) = 120
#define B2_MAX_BORDERS ( B2_CLUSTER_COUNT * ( B2_CLUSTER_COUNT - 1 ) / 2 )

typedef struct b2Cluster
{
	b2IntArray bodyIndices;
	b2Vec2 center;
	b2Vec2 accumulator;
} b2Cluster;

// Interior constraints for one cluster, allocated from the arena each step
typedef struct b2ClusterSolveData
{
	b2ContactSim** contacts;
	int contactCount;

	b2JointSim** joints;
	int jointCount;

	b2ContactConstraint* contactConstraints;

	// Signaled by the worker when this cluster's solve phase is done
	b2AtomicInt solveComplete;
} b2ClusterSolveData;

// Border constraints between two clusters (clusterA < clusterB)
typedef struct b2BorderConstraints
{
	int clusterA;
	int clusterB;

	b2ContactSim** contacts;
	int contactCount;

	b2JointSim** joints;
	int jointCount;

	b2ContactConstraint* contactConstraints;
} b2BorderConstraints;

typedef struct b2ClusterManager
{
	b2Cluster clusters[B2_CLUSTER_COUNT];
} b2ClusterManager;

void b2CreateClusters( b2ClusterManager* manager );
void b2DestroyClusters( b2ClusterManager* manager );

void b2ComputeClusters( b2World* world );

// Classify all awake touching contacts and joints into cluster interiors and borders.
// Allocates from the arena; data is transient for this step.
void b2ClassifyConstraints( b2World* world, b2StepContext* context );
