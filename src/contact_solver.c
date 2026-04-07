// SPDX-FileCopyrightText: 2023 Erin Catto
// SPDX-License-Identifier: MIT

#include "contact_solver.h"

#include "body.h"
#include "contact.h"
#include "core.h"
#include "physics_world.h"

#include <stddef.h>

void b2PrepareContactConstraints( b2StepContext* context, b2ContactSim** contacts, b2ContactConstraint* constraints, int count )
{
	b2World* world = context->world;
	b2Body* bodies = world->bodies.data;
	float warmStartScale = world->enableWarmStarting ? 1.0f : 0.0f;

	for ( int i = 0; i < count; ++i )
	{
		b2ContactSim* contactSim = contacts[i];
		b2Body* bodyA = bodies + contactSim->bodyIdA;
		b2Body* bodyB = bodies + contactSim->bodyIdB;
		int indexA = bodyA->stateIndex;
		int indexB = bodyB->stateIndex;

		const b2Manifold* manifold = &contactSim->manifold;
		int pointCount = manifold->pointCount;
		B2_ASSERT( 0 < pointCount && pointCount <= 2 );

		b2ContactConstraint* constraint = constraints + i;
		constraint->indexA = indexA;
		constraint->indexB = indexB;
		constraint->normal = manifold->normal;
		constraint->friction = contactSim->friction;
		constraint->restitution = contactSim->restitution;
		constraint->rollingResistance = contactSim->rollingResistance;
		constraint->rollingImpulse = warmStartScale * manifold->rollingImpulse;
		constraint->tangentSpeed = contactSim->tangentSpeed;
		constraint->pointCount = pointCount;

		b2Vec2 vA = bodyA->linearVelocity;
		float wA = bodyA->angularVelocity;
		float mA = bodyA->invMass;
		float iA = bodyA->invInertia;

		b2Vec2 vB = bodyB->linearVelocity;
		float wB = bodyB->angularVelocity;
		float mB = bodyB->invMass;
		float iB = bodyB->invInertia;

		{
			float k = iA + iB;
			constraint->rollingMass = k > 0.0f ? 1.0f / k : 0.0f;
		}

		b2Vec2 normal = constraint->normal;
		b2Vec2 tangent = b2RightPerp( normal );

		for ( int j = 0; j < pointCount; ++j )
		{
			const b2ManifoldPoint* mp = manifold->points + j;
			b2ContactConstraintPoint* cp = constraint->points + j;

			cp->normalImpulse = warmStartScale * mp->normalImpulse;
			cp->tangentImpulse = warmStartScale * mp->tangentImpulse;
			cp->totalNormalImpulse = 0.0f;

			b2Vec2 rA = mp->anchorA;
			b2Vec2 rB = mp->anchorB;
			cp->anchorA = rA;
			cp->anchorB = rB;
			cp->baseSeparation = mp->separation - b2Dot( b2Sub( rB, rA ), normal );

			float rnA = b2Cross( rA, normal );
			float rnB = b2Cross( rB, normal );
			float kNormal = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			cp->normalMass = kNormal > 0.0f ? 1.0f / kNormal : 0.0f;

			float rtA = b2Cross( rA, tangent );
			float rtB = b2Cross( rB, tangent );
			float kTangent = mA + mB + iA * rtA * rtA + iB * rtB * rtB;
			cp->tangentMass = kTangent > 0.0f ? 1.0f / kTangent : 0.0f;

			b2Vec2 vrA = b2Add( vA, b2CrossSV( wA, rA ) );
			b2Vec2 vrB = b2Add( vB, b2CrossSV( wB, rB ) );
			cp->relativeVelocity = b2Dot( normal, b2Sub( vrB, vrA ) );
		}
	}
}

void b2WarmStartContactConstraints( b2StepContext* context, b2ContactConstraint* constraints, int count )
{
	b2BodyState* states = context->states;
	b2BodyState dummyState = b2_identityBodyState;

	for ( int i = 0; i < count; ++i )
	{
		b2ContactConstraint* constraint = constraints + i;

		int indexA = constraint->indexA;
		int indexB = constraint->indexB;

		b2BodyState* stateA = indexA == B2_NULL_INDEX ? &dummyState : states + indexA;
		b2BodyState* stateB = indexB == B2_NULL_INDEX ? &dummyState : states + indexB;

		b2Vec2 vA = stateA->linearVelocity;
		float wA = stateA->angularVelocity;
		float mA = stateA->invMass;
		float iA = stateA->invInertia;

		b2Vec2 vB = stateB->linearVelocity;
		float wB = stateB->angularVelocity;
		float mB = stateB->invMass;
		float iB = stateB->invInertia;

		b2Vec2 normal = constraint->normal;
		b2Vec2 tangent = b2RightPerp( normal );
		int pointCount = constraint->pointCount;

		for ( int j = 0; j < pointCount; ++j )
		{
			b2ContactConstraintPoint* cp = constraint->points + j;
			b2Vec2 rA = cp->anchorA;
			b2Vec2 rB = cp->anchorB;

			b2Vec2 P = b2Add( b2MulSV( cp->normalImpulse, normal ), b2MulSV( cp->tangentImpulse, tangent ) );
			cp->totalNormalImpulse += cp->normalImpulse;

			wA -= iA * b2Cross( rA, P );
			vA = b2MulAdd( vA, -mA, P );
			wB += iB * b2Cross( rB, P );
			vB = b2MulAdd( vB, mB, P );
		}

		wA -= iA * constraint->rollingImpulse;
		wB += iB * constraint->rollingImpulse;

		if ( stateA->flags & b2_dynamicFlag )
		{
			stateA->linearVelocity = vA;
			stateA->angularVelocity = wA;
		}

		if ( stateB->flags & b2_dynamicFlag )
		{
			stateB->linearVelocity = vB;
			stateB->angularVelocity = wB;
		}
	}
}

void b2SolveContactConstraints( b2StepContext* context, b2ContactConstraint* constraints, int count, float inv_h,
								float contactSpeed, bool useBias )
{
	b2BodyState dummyState = b2_identityBodyState;

	b2Softness contactSoftness = context->contactSoftness;
	b2Softness staticSoftness = context->staticSoftness;
	b2BodyState* states = context->states;

	for ( int i = 0; i < count; ++i )
	{
		b2ContactConstraint* constraint = constraints + i;

		int indexA = constraint->indexA;
		int indexB = constraint->indexB;

		b2BodyState* stateA = indexA == B2_NULL_INDEX ? &dummyState : states + indexA;
		b2Vec2 vA = stateA->linearVelocity;
		float wA = stateA->angularVelocity;
		float mA = stateA->invMass;
		float iA = stateA->invInertia;
		b2Rot dqA = stateA->deltaRotation;

		b2BodyState* stateB = indexB == B2_NULL_INDEX ? &dummyState : states + indexB;
		b2Vec2 vB = stateB->linearVelocity;
		float wB = stateB->angularVelocity;
		float mB = stateB->invMass;
		float iB = stateB->invInertia;
		b2Rot dqB = stateB->deltaRotation;

		b2Vec2 dp = b2Sub( stateB->deltaPosition, stateA->deltaPosition );

		b2Vec2 normal = constraint->normal;
		b2Vec2 tangent = b2RightPerp( normal );
		float friction = constraint->friction;
		b2Softness softness = ( indexA == B2_NULL_INDEX || indexB == B2_NULL_INDEX ) ? staticSoftness : contactSoftness;

		int pointCount = constraint->pointCount;
		float totalNormalImpulse = 0.0f;

		// Non-penetration
		for ( int j = 0; j < pointCount; ++j )
		{
			b2ContactConstraintPoint* cp = constraint->points + j;
			b2Vec2 rA = cp->anchorA;
			b2Vec2 rB = cp->anchorB;

			b2Vec2 ds = b2Add( dp, b2Sub( b2RotateVector( dqB, rB ), b2RotateVector( dqA, rA ) ) );
			float s = cp->baseSeparation + b2Dot( ds, normal );

			float velocityBias = 0.0f;
			float massScale = 1.0f;
			float impulseScale = 0.0f;
			if ( s > 0.0f )
			{
				velocityBias = s * inv_h;
			}
			else if ( useBias )
			{
				velocityBias = b2MaxFloat( softness.massScale * softness.biasRate * s, -contactSpeed );
				massScale = softness.massScale;
				impulseScale = softness.impulseScale;
			}

			b2Vec2 vrA = b2Add( vA, b2CrossSV( wA, rA ) );
			b2Vec2 vrB = b2Add( vB, b2CrossSV( wB, rB ) );
			float vn = b2Dot( b2Sub( vrB, vrA ), normal );

			float impulse = -cp->normalMass * ( massScale * vn + velocityBias ) - impulseScale * cp->normalImpulse;
			float newImpulse = b2MaxFloat( cp->normalImpulse + impulse, 0.0f );
			impulse = newImpulse - cp->normalImpulse;
			cp->normalImpulse = newImpulse;
			cp->totalNormalImpulse += impulse;

			totalNormalImpulse += newImpulse;

			b2Vec2 P = b2MulSV( impulse, normal );
			vA = b2MulSub( vA, mA, P );
			wA -= iA * b2Cross( rA, P );
			vB = b2MulAdd( vB, mB, P );
			wB += iB * b2Cross( rB, P );
		}

		// Friction
		for ( int j = 0; j < pointCount; ++j )
		{
			b2ContactConstraintPoint* cp = constraint->points + j;
			b2Vec2 rA = cp->anchorA;
			b2Vec2 rB = cp->anchorB;

			b2Vec2 vrB = b2Add( vB, b2CrossSV( wB, rB ) );
			b2Vec2 vrA = b2Add( vA, b2CrossSV( wA, rA ) );
			float vt = b2Dot( b2Sub( vrB, vrA ), tangent ) - constraint->tangentSpeed;

			float impulse = cp->tangentMass * ( -vt );
			float maxFriction = friction * cp->normalImpulse;
			float newImpulse = b2ClampFloat( cp->tangentImpulse + impulse, -maxFriction, maxFriction );
			impulse = newImpulse - cp->tangentImpulse;
			cp->tangentImpulse = newImpulse;

			b2Vec2 P = b2MulSV( impulse, tangent );
			vA = b2MulSub( vA, mA, P );
			wA -= iA * b2Cross( rA, P );
			vB = b2MulAdd( vB, mB, P );
			wB += iB * b2Cross( rB, P );
		}

		// Rolling resistance
		if ( constraint->rollingResistance > 0.0f )
		{
			float deltaLambda = -constraint->rollingMass * ( wB - wA );
			float lambda = constraint->rollingImpulse;
			float maxLambda = constraint->rollingResistance * totalNormalImpulse;
			constraint->rollingImpulse = b2ClampFloat( lambda + deltaLambda, -maxLambda, maxLambda );
			deltaLambda = constraint->rollingImpulse - lambda;

			wA -= iA * deltaLambda;
			wB += iB * deltaLambda;
		}

		if ( stateA->flags & b2_dynamicFlag )
		{
			stateA->linearVelocity = vA;
			stateA->angularVelocity = wA;
		}

		if ( stateB->flags & b2_dynamicFlag )
		{
			stateB->linearVelocity = vB;
			stateB->angularVelocity = wB;
		}
	}
}

void b2ApplyContactRestitution( b2StepContext* context, b2ContactConstraint* constraints, int count, float threshold )
{
	b2BodyState* states = context->states;
	b2BodyState dummyState = b2_identityBodyState;

	for ( int i = 0; i < count; ++i )
	{
		b2ContactConstraint* constraint = constraints + i;
		float restitution = constraint->restitution;
		if ( restitution == 0.0f )
		{
			continue;
		}

		int indexA = constraint->indexA;
		int indexB = constraint->indexB;

		b2BodyState* stateA = indexA == B2_NULL_INDEX ? &dummyState : states + indexA;
		b2Vec2 vA = stateA->linearVelocity;
		float wA = stateA->angularVelocity;
		float mA = stateA->invMass;
		float iA = stateA->invInertia;

		b2BodyState* stateB = indexB == B2_NULL_INDEX ? &dummyState : states + indexB;
		b2Vec2 vB = stateB->linearVelocity;
		float wB = stateB->angularVelocity;
		float mB = stateB->invMass;
		float iB = stateB->invInertia;

		b2Vec2 normal = constraint->normal;
		int pointCount = constraint->pointCount;

		for ( int j = 0; j < pointCount; ++j )
		{
			b2ContactConstraintPoint* cp = constraint->points + j;

			if ( cp->relativeVelocity > -threshold || cp->totalNormalImpulse == 0.0f )
			{
				continue;
			}

			b2Vec2 rA = cp->anchorA;
			b2Vec2 rB = cp->anchorB;

			b2Vec2 vrB = b2Add( vB, b2CrossSV( wB, rB ) );
			b2Vec2 vrA = b2Add( vA, b2CrossSV( wA, rA ) );
			float vn = b2Dot( b2Sub( vrB, vrA ), normal );

			float impulse = -cp->normalMass * ( vn + restitution * cp->relativeVelocity );
			float newImpulse = b2MaxFloat( cp->normalImpulse + impulse, 0.0f );
			impulse = newImpulse - cp->normalImpulse;
			cp->normalImpulse = newImpulse;
			cp->totalNormalImpulse += impulse;

			b2Vec2 P = b2MulSV( impulse, normal );
			vA = b2MulSub( vA, mA, P );
			wA -= iA * b2Cross( rA, P );
			vB = b2MulAdd( vB, mB, P );
			wB += iB * b2Cross( rB, P );
		}

		if ( stateA->flags & b2_dynamicFlag )
		{
			stateA->linearVelocity = vA;
			stateA->angularVelocity = wA;
		}

		if ( stateB->flags & b2_dynamicFlag )
		{
			stateB->linearVelocity = vB;
			stateB->angularVelocity = wB;
		}
	}
}

void b2StoreContactImpulses( b2ContactSim** contacts, b2ContactConstraint* constraints, int count )
{
	for ( int i = 0; i < count; ++i )
	{
		const b2ContactConstraint* constraint = constraints + i;
		b2ContactSim* contact = contacts[i];
		b2Manifold* manifold = &contact->manifold;
		int pointCount = manifold->pointCount;

		for ( int j = 0; j < pointCount; ++j )
		{
			manifold->points[j].normalImpulse = constraint->points[j].normalImpulse;
			manifold->points[j].tangentImpulse = constraint->points[j].tangentImpulse;
			manifold->points[j].totalNormalImpulse = constraint->points[j].totalNormalImpulse;
			manifold->points[j].normalVelocity = constraint->points[j].relativeVelocity;
		}

		manifold->rollingImpulse = constraint->rollingImpulse;
	}
}

#if 0
#if defined( B2_SIMD_AVX2 )

#include <immintrin.h>

// wide float holds 8 numbers
typedef __m256 b2FloatW;

#elif defined( B2_SIMD_NEON )

#include <arm_neon.h>

// wide float holds 4 numbers
typedef float32x4_t b2FloatW;

#elif defined( B2_SIMD_SSE2 )

#include <emmintrin.h>

// wide float holds 4 numbers
typedef __m128 b2FloatW;

#else

// scalar math
typedef struct b2FloatW
{
	float x, y, z, w;
} b2FloatW;

#endif

// Wide vec2
typedef struct b2Vec2W
{
	b2FloatW X, Y;
} b2Vec2W;

// Wide rotation
typedef struct b2RotW
{
	b2FloatW C, S;
} b2RotW;

#if defined( B2_SIMD_AVX2 )

static inline b2FloatW b2ZeroW( void )
{
	return _mm256_setzero_ps();
}

static inline b2FloatW b2SplatW( float scalar )
{
	return _mm256_set1_ps( scalar );
}

static inline b2FloatW b2AddW( b2FloatW a, b2FloatW b )
{
	return _mm256_add_ps( a, b );
}

static inline b2FloatW b2SubW( b2FloatW a, b2FloatW b )
{
	return _mm256_sub_ps( a, b );
}

static inline b2FloatW b2MulW( b2FloatW a, b2FloatW b )
{
	return _mm256_mul_ps( a, b );
}

static inline b2FloatW b2MulAddW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	// FMA can be emulated: https://github.com/lattera/glibc/blob/master/sysdeps/ieee754/dbl-64/s_fmaf.c#L34
	// return _mm256_fmadd_ps( b, c, a );
	return _mm256_add_ps( _mm256_mul_ps( b, c ), a );
}

static inline b2FloatW b2MulSubW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	// return _mm256_fnmadd_ps(b, c, a);
	return _mm256_sub_ps( a, _mm256_mul_ps( b, c ) );
}

static inline b2FloatW b2MinW( b2FloatW a, b2FloatW b )
{
	return _mm256_min_ps( a, b );
}

static inline b2FloatW b2MaxW( b2FloatW a, b2FloatW b )
{
	return _mm256_max_ps( a, b );
}

// a = clamp(a, -b, b)
static inline b2FloatW b2SymClampW( b2FloatW a, b2FloatW b )
{
	b2FloatW nb = _mm256_sub_ps( _mm256_setzero_ps(), b );
	return _mm256_max_ps( nb, _mm256_min_ps( a, b ) );
}

static inline b2FloatW b2OrW( b2FloatW a, b2FloatW b )
{
	return _mm256_or_ps( a, b );
}

static inline b2FloatW b2GreaterThanW( b2FloatW a, b2FloatW b )
{
	return _mm256_cmp_ps( a, b, _CMP_GT_OQ );
}

static inline b2FloatW b2EqualsW( b2FloatW a, b2FloatW b )
{
	return _mm256_cmp_ps( a, b, _CMP_EQ_OQ );
}

static inline bool b2AllZeroW( b2FloatW a )
{
	// Compare each element with zero
	b2FloatW zero = _mm256_setzero_ps();
	b2FloatW cmp = _mm256_cmp_ps( a, zero, _CMP_EQ_OQ );

	// Create a mask from the comparison results
	int mask = _mm256_movemask_ps( cmp );

	// If all elements are zero, the mask will be 0xFF (11111111 in binary)
	return mask == 0xFF;
}

// component-wise returns mask ? b : a
static inline b2FloatW b2BlendW( b2FloatW a, b2FloatW b, b2FloatW mask )
{
	return _mm256_blendv_ps( a, b, mask );
}

#elif defined( B2_SIMD_NEON )

static inline b2FloatW b2ZeroW( void )
{
	return vdupq_n_f32( 0.0f );
}

static inline b2FloatW b2SplatW( float scalar )
{
	return vdupq_n_f32( scalar );
}

static inline b2FloatW b2SetW( float a, float b, float c, float d )
{
	float32_t array[4] = { a, b, c, d };
	return vld1q_f32( array );
}

static inline b2FloatW b2AddW( b2FloatW a, b2FloatW b )
{
	return vaddq_f32( a, b );
}

static inline b2FloatW b2SubW( b2FloatW a, b2FloatW b )
{
	return vsubq_f32( a, b );
}

static inline b2FloatW b2MulW( b2FloatW a, b2FloatW b )
{
	return vmulq_f32( a, b );
}

static inline b2FloatW b2MulAddW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	return vaddq_f32( a, vmulq_f32( b, c ) );
}

static inline b2FloatW b2MulSubW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	return vsubq_f32( a, vmulq_f32( b, c ) );
}

static inline b2FloatW b2MinW( b2FloatW a, b2FloatW b )
{
	return vminq_f32( a, b );
}

static inline b2FloatW b2MaxW( b2FloatW a, b2FloatW b )
{
	return vmaxq_f32( a, b );
}

// a = clamp(a, -b, b)
static inline b2FloatW b2SymClampW( b2FloatW a, b2FloatW b )
{
	b2FloatW nb = vnegq_f32( b );
	return vmaxq_f32( nb, vminq_f32( a, b ) );
}

static inline b2FloatW b2OrW( b2FloatW a, b2FloatW b )
{
	return vreinterpretq_f32_u32( vorrq_u32( vreinterpretq_u32_f32( a ), vreinterpretq_u32_f32( b ) ) );
}

static inline b2FloatW b2GreaterThanW( b2FloatW a, b2FloatW b )
{
	return vreinterpretq_f32_u32( vcgtq_f32( a, b ) );
}

static inline b2FloatW b2EqualsW( b2FloatW a, b2FloatW b )
{
	return vreinterpretq_f32_u32( vceqq_f32( a, b ) );
}

static inline bool b2AllZeroW( b2FloatW a )
{
	// Create a zero vector for comparison
	b2FloatW zero = vdupq_n_f32( 0.0f );

	// Compare the input vector with zero
	uint32x4_t cmp_result = vceqq_f32( a, zero );

// Check if all comparison results are non-zero using vminvq
#ifdef __ARM_FEATURE_SVE
	// ARM v8.2+ has horizontal minimum instruction
	return vminvq_u32( cmp_result ) != 0;
#else
	// For older ARM architectures, we need to manually check all lanes
	return vgetq_lane_u32( cmp_result, 0 ) != 0 && vgetq_lane_u32( cmp_result, 1 ) != 0 && vgetq_lane_u32( cmp_result, 2 ) != 0 &&
		   vgetq_lane_u32( cmp_result, 3 ) != 0;
#endif
}

// component-wise returns mask ? b : a
static inline b2FloatW b2BlendW( b2FloatW a, b2FloatW b, b2FloatW mask )
{
	uint32x4_t mask32 = vreinterpretq_u32_f32( mask );
	return vbslq_f32( mask32, b, a );
}

static inline b2FloatW b2LoadW( const float32_t* data )
{
	return vld1q_f32( data );
}

static inline void b2StoreW( float32_t* data, b2FloatW a )
{
	vst1q_f32( data, a );
}

static inline b2FloatW b2UnpackLoW( b2FloatW a, b2FloatW b )
{
#if defined( _M_ARM64 ) || defined( __aarch64__ )
	return vzip1q_f32( a, b );
#else
	float32x2_t a1 = vget_low_f32( a );
	float32x2_t b1 = vget_low_f32( b );
	float32x2x2_t result = vzip_f32( a1, b1 );
	return vcombine_f32( result.val[0], result.val[1] );
#endif
}

static inline b2FloatW b2UnpackHiW( b2FloatW a, b2FloatW b )
{
#if defined( _M_ARM64 ) || defined( __aarch64__ )
	return vzip2q_f32( a, b );
#else
	float32x2_t a1 = vget_high_f32( a );
	float32x2_t b1 = vget_high_f32( b );
	float32x2x2_t result = vzip_f32( a1, b1 );
	return vcombine_f32( result.val[0], result.val[1] );
#endif
}

#elif defined( B2_SIMD_SSE2 )

static inline b2FloatW b2ZeroW( void )
{
	return _mm_setzero_ps();
}

static inline b2FloatW b2SplatW( float scalar )
{
	return _mm_set1_ps( scalar );
}

static inline b2FloatW b2SetW( float a, float b, float c, float d )
{
	return _mm_setr_ps( a, b, c, d );
}

static inline b2FloatW b2AddW( b2FloatW a, b2FloatW b )
{
	return _mm_add_ps( a, b );
}

static inline b2FloatW b2SubW( b2FloatW a, b2FloatW b )
{
	return _mm_sub_ps( a, b );
}

static inline b2FloatW b2MulW( b2FloatW a, b2FloatW b )
{
	return _mm_mul_ps( a, b );
}

static inline b2FloatW b2MulAddW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	return _mm_add_ps( a, _mm_mul_ps( b, c ) );
}

static inline b2FloatW b2MulSubW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	return _mm_sub_ps( a, _mm_mul_ps( b, c ) );
}

static inline b2FloatW b2MinW( b2FloatW a, b2FloatW b )
{
	return _mm_min_ps( a, b );
}

static inline b2FloatW b2MaxW( b2FloatW a, b2FloatW b )
{
	return _mm_max_ps( a, b );
}

// a = clamp(a, -b, b)
static inline b2FloatW b2SymClampW( b2FloatW a, b2FloatW b )
{
	// Create a mask with the sign bit set for each element
	__m128 mask = _mm_set1_ps( -0.0f );

	// XOR the input with the mask to negate each element
	__m128 nb = _mm_xor_ps( b, mask );

	return _mm_max_ps( nb, _mm_min_ps( a, b ) );
}

static inline b2FloatW b2OrW( b2FloatW a, b2FloatW b )
{
	return _mm_or_ps( a, b );
}

static inline b2FloatW b2GreaterThanW( b2FloatW a, b2FloatW b )
{
	return _mm_cmpgt_ps( a, b );
}

static inline b2FloatW b2EqualsW( b2FloatW a, b2FloatW b )
{
	return _mm_cmpeq_ps( a, b );
}

static inline bool b2AllZeroW( b2FloatW a )
{
	// Compare each element with zero
	b2FloatW zero = _mm_setzero_ps();
	b2FloatW cmp = _mm_cmpeq_ps( a, zero );

	// Create a mask from the comparison results
	int mask = _mm_movemask_ps( cmp );

	// If all elements are zero, the mask will be 0xF (1111 in binary)
	return mask == 0xF;
}

// component-wise returns mask ? b : a
static inline b2FloatW b2BlendW( b2FloatW a, b2FloatW b, b2FloatW mask )
{
	return _mm_or_ps( _mm_and_ps( mask, b ), _mm_andnot_ps( mask, a ) );
}

static inline b2FloatW b2LoadW( const float* data )
{
	return _mm_load_ps( data );
}

static inline void b2StoreW( float* data, b2FloatW a )
{
	_mm_store_ps( data, a );
}

static inline b2FloatW b2UnpackLoW( b2FloatW a, b2FloatW b )
{
	return _mm_unpacklo_ps( a, b );
}

static inline b2FloatW b2UnpackHiW( b2FloatW a, b2FloatW b )
{
	return _mm_unpackhi_ps( a, b );
}

#else

static inline b2FloatW b2ZeroW( void )
{
	return (b2FloatW){ 0.0f, 0.0f, 0.0f, 0.0f };
}

static inline b2FloatW b2SplatW( float scalar )
{
	return (b2FloatW){ scalar, scalar, scalar, scalar };
}

static inline b2FloatW b2AddW( b2FloatW a, b2FloatW b )
{
	return (b2FloatW){ a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w };
}

static inline b2FloatW b2SubW( b2FloatW a, b2FloatW b )
{
	return (b2FloatW){ a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w };
}

static inline b2FloatW b2MulW( b2FloatW a, b2FloatW b )
{
	return (b2FloatW){ a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w };
}

static inline b2FloatW b2MulAddW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	return (b2FloatW){ a.x + b.x * c.x, a.y + b.y * c.y, a.z + b.z * c.z, a.w + b.w * c.w };
}

static inline b2FloatW b2MulSubW( b2FloatW a, b2FloatW b, b2FloatW c )
{
	return (b2FloatW){ a.x - b.x * c.x, a.y - b.y * c.y, a.z - b.z * c.z, a.w - b.w * c.w };
}

static inline b2FloatW b2MinW( b2FloatW a, b2FloatW b )
{
	b2FloatW r;
	r.x = a.x <= b.x ? a.x : b.x;
	r.y = a.y <= b.y ? a.y : b.y;
	r.z = a.z <= b.z ? a.z : b.z;
	r.w = a.w <= b.w ? a.w : b.w;
	return r;
}

static inline b2FloatW b2MaxW( b2FloatW a, b2FloatW b )
{
	b2FloatW r;
	r.x = a.x >= b.x ? a.x : b.x;
	r.y = a.y >= b.y ? a.y : b.y;
	r.z = a.z >= b.z ? a.z : b.z;
	r.w = a.w >= b.w ? a.w : b.w;
	return r;
}

// a = clamp(a, -b, b)
static inline b2FloatW b2SymClampW( b2FloatW a, b2FloatW b )
{
	b2FloatW r;
	r.x = b2ClampFloat( a.x, -b.x, b.x );
	r.y = b2ClampFloat( a.y, -b.y, b.y );
	r.z = b2ClampFloat( a.z, -b.z, b.z );
	r.w = b2ClampFloat( a.w, -b.w, b.w );
	return r;
}

static inline b2FloatW b2OrW( b2FloatW a, b2FloatW b )
{
	b2FloatW r;
	r.x = a.x != 0.0f || b.x != 0.0f ? 1.0f : 0.0f;
	r.y = a.y != 0.0f || b.y != 0.0f ? 1.0f : 0.0f;
	r.z = a.z != 0.0f || b.z != 0.0f ? 1.0f : 0.0f;
	r.w = a.w != 0.0f || b.w != 0.0f ? 1.0f : 0.0f;
	return r;
}

static inline b2FloatW b2GreaterThanW( b2FloatW a, b2FloatW b )
{
	b2FloatW r;
	r.x = a.x > b.x ? 1.0f : 0.0f;
	r.y = a.y > b.y ? 1.0f : 0.0f;
	r.z = a.z > b.z ? 1.0f : 0.0f;
	r.w = a.w > b.w ? 1.0f : 0.0f;
	return r;
}

static inline b2FloatW b2EqualsW( b2FloatW a, b2FloatW b )
{
	b2FloatW r;
	r.x = a.x == b.x ? 1.0f : 0.0f;
	r.y = a.y == b.y ? 1.0f : 0.0f;
	r.z = a.z == b.z ? 1.0f : 0.0f;
	r.w = a.w == b.w ? 1.0f : 0.0f;
	return r;
}

static inline bool b2AllZeroW( b2FloatW a )
{
	return a.x == 0.0f && a.y == 0.0f && a.z == 0.0f && a.w == 0.0f;
}

// component-wise returns mask ? b : a
static inline b2FloatW b2BlendW( b2FloatW a, b2FloatW b, b2FloatW mask )
{
	b2FloatW r;
	r.x = mask.x != 0.0f ? b.x : a.x;
	r.y = mask.y != 0.0f ? b.y : a.y;
	r.z = mask.z != 0.0f ? b.z : a.z;
	r.w = mask.w != 0.0f ? b.w : a.w;
	return r;
}

#endif

static inline b2FloatW b2DotW( b2Vec2W a, b2Vec2W b )
{
	return b2AddW( b2MulW( a.X, b.X ), b2MulW( a.Y, b.Y ) );
}

static inline b2FloatW b2CrossW( b2Vec2W a, b2Vec2W b )
{
	return b2SubW( b2MulW( a.X, b.Y ), b2MulW( a.Y, b.X ) );
}

static inline b2Vec2W b2RotateVectorW( b2RotW q, b2Vec2W v )
{
	return (b2Vec2W){ b2SubW( b2MulW( q.C, v.X ), b2MulW( q.S, v.Y ) ), b2AddW( b2MulW( q.S, v.X ), b2MulW( q.C, v.Y ) ) };
}

// Soft contact constraints with sub-stepping support
// Uses fixed anchors for Jacobians for better behavior on rolling shapes (circles & capsules)
// http://mmacklin.com/smallsteps.pdf
// https://box2d.org/files/ErinCatto_SoftConstraints_GDC2011.pdf

typedef struct b2ContactConstraintWide
{
	int indexA[B2_SIMD_WIDTH];
	int indexB[B2_SIMD_WIDTH];

	b2FloatW invMassA, invMassB;
	b2FloatW invIA, invIB;
	b2Vec2W normal;
	b2FloatW friction;
	b2FloatW tangentSpeed;
	b2FloatW rollingResistance;
	b2FloatW rollingMass;
	b2FloatW rollingImpulse;
	b2FloatW biasRate;
	b2FloatW massScale;
	b2FloatW impulseScale;
	b2Vec2W anchorA1, anchorB1;
	b2FloatW normalMass1, tangentMass1;
	b2FloatW baseSeparation1;
	b2FloatW normalImpulse1;
	b2FloatW totalNormalImpulse1;
	b2FloatW tangentImpulse1;
	b2Vec2W anchorA2, anchorB2;
	b2FloatW baseSeparation2;
	b2FloatW normalImpulse2;
	b2FloatW totalNormalImpulse2;
	b2FloatW tangentImpulse2;
	b2FloatW normalMass2, tangentMass2;
	b2FloatW restitution;
	b2FloatW relativeVelocity1, relativeVelocity2;
} b2ContactConstraintWide;

int b2GetContactConstraintSIMDByteCount( void )
{
	return sizeof( b2ContactConstraintWide );
}

// wide version of b2BodyState
typedef struct b2BodyStateW
{
	b2Vec2W v;
	b2FloatW w;
	b2FloatW flags;
	b2Vec2W dp;
	b2RotW dq;
} b2BodyStateW;

// Custom gather/scatter for each SIMD type
#if defined( B2_SIMD_AVX2 )

// This is a load and 8x8 transpose
static b2BodyStateW b2GatherBodies( const b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices )
{
	_Static_assert( sizeof( b2BodyState ) == 32, "b2BodyState not 32 bytes" );
	B2_ASSERT( ( (uintptr_t)states & 0x1F ) == 0 );

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;
	int i5 = indices[4] - 1;
	int i6 = indices[5] - 1;
	int i7 = indices[6] - 1;
	int i8 = indices[7] - 1;

	// b2BodyState b2_identityBodyState = {{0.0f, 0.0f}, 0.0f, 0, {0.0f, 0.0f}, {1.0f, 0.0f}};
	b2FloatW identity = _mm256_setr_ps( 0.0f, 0.0f, 0.0f, 0, 0.0f, 0.0f, 1.0f, 0.0f );
	b2FloatW b0 = i1 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i1 ) );
	b2FloatW b1 = i2 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i2 ) );
	b2FloatW b2 = i3 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i3 ) );
	b2FloatW b3 = i4 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i4 ) );
	b2FloatW b4 = i5 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i5 ) );
	b2FloatW b5 = i6 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i6 ) );
	b2FloatW b6 = i7 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i7 ) );
	b2FloatW b7 = i8 == B2_NULL_INDEX ? identity : _mm256_load_ps( (float*)( states + i8 ) );

	b2FloatW t0 = _mm256_unpacklo_ps( b0, b1 );
	b2FloatW t1 = _mm256_unpackhi_ps( b0, b1 );
	b2FloatW t2 = _mm256_unpacklo_ps( b2, b3 );
	b2FloatW t3 = _mm256_unpackhi_ps( b2, b3 );
	b2FloatW t4 = _mm256_unpacklo_ps( b4, b5 );
	b2FloatW t5 = _mm256_unpackhi_ps( b4, b5 );
	b2FloatW t6 = _mm256_unpacklo_ps( b6, b7 );
	b2FloatW t7 = _mm256_unpackhi_ps( b6, b7 );
	b2FloatW tt0 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt1 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 3, 2, 3, 2 ) );
	b2FloatW tt2 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt3 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );
	b2FloatW tt4 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt5 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 3, 2, 3, 2 ) );
	b2FloatW tt6 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt7 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 3, 2, 3, 2 ) );

	b2BodyStateW simdBody;
	simdBody.v.X = _mm256_permute2f128_ps( tt0, tt4, 0x20 );
	simdBody.v.Y = _mm256_permute2f128_ps( tt1, tt5, 0x20 );
	simdBody.w = _mm256_permute2f128_ps( tt2, tt6, 0x20 );
	simdBody.flags = _mm256_permute2f128_ps( tt3, tt7, 0x20 );
	simdBody.dp.X = _mm256_permute2f128_ps( tt0, tt4, 0x31 );
	simdBody.dp.Y = _mm256_permute2f128_ps( tt1, tt5, 0x31 );
	simdBody.dq.C = _mm256_permute2f128_ps( tt2, tt6, 0x31 );
	simdBody.dq.S = _mm256_permute2f128_ps( tt3, tt7, 0x31 );
	return simdBody;
}

// This writes everything back to the solver bodies but only the velocities change
static void b2ScatterBodies( b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices, const b2BodyStateW* B2_RESTRICT simdBody )
{
	_Static_assert( sizeof( b2BodyState ) == 32, "b2BodyState not 32 bytes" );
	B2_ASSERT( ( (uintptr_t)states & 0x1F ) == 0 );
	b2FloatW t0 = _mm256_unpacklo_ps( simdBody->v.X, simdBody->v.Y );
	b2FloatW t1 = _mm256_unpackhi_ps( simdBody->v.X, simdBody->v.Y );
	b2FloatW t2 = _mm256_unpacklo_ps( simdBody->w, simdBody->flags );
	b2FloatW t3 = _mm256_unpackhi_ps( simdBody->w, simdBody->flags );
	b2FloatW t4 = _mm256_unpacklo_ps( simdBody->dp.X, simdBody->dp.Y );
	b2FloatW t5 = _mm256_unpackhi_ps( simdBody->dp.X, simdBody->dp.Y );
	b2FloatW t6 = _mm256_unpacklo_ps( simdBody->dq.C, simdBody->dq.S );
	b2FloatW t7 = _mm256_unpackhi_ps( simdBody->dq.C, simdBody->dq.S );
	b2FloatW tt0 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt1 = _mm256_shuffle_ps( t0, t2, _MM_SHUFFLE( 3, 2, 3, 2 ) );
	b2FloatW tt2 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt3 = _mm256_shuffle_ps( t1, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) );
	b2FloatW tt4 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt5 = _mm256_shuffle_ps( t4, t6, _MM_SHUFFLE( 3, 2, 3, 2 ) );
	b2FloatW tt6 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 1, 0, 1, 0 ) );
	b2FloatW tt7 = _mm256_shuffle_ps( t5, t7, _MM_SHUFFLE( 3, 2, 3, 2 ) );

	// I don't use any dummy body in the body array because this will lead to multithreaded sharing and the
	// associated cache flushing.

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;
	int i5 = indices[4] - 1;
	int i6 = indices[5] - 1;
	int i7 = indices[6] - 1;
	int i8 = indices[7] - 1;

	if ( i1 != B2_NULL_INDEX && ( states[i1].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i1 ), _mm256_permute2f128_ps( tt0, tt4, 0x20 ) );
	if ( i2 != B2_NULL_INDEX && ( states[i2].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i2 ), _mm256_permute2f128_ps( tt1, tt5, 0x20 ) );
	if ( i3 != B2_NULL_INDEX && ( states[i3].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i3 ), _mm256_permute2f128_ps( tt2, tt6, 0x20 ) );
	if ( i4 != B2_NULL_INDEX && ( states[i4].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i4 ), _mm256_permute2f128_ps( tt3, tt7, 0x20 ) );
	if ( i5 != B2_NULL_INDEX && ( states[i5].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i5 ), _mm256_permute2f128_ps( tt0, tt4, 0x31 ) );
	if ( i6 != B2_NULL_INDEX && ( states[i6].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i6 ), _mm256_permute2f128_ps( tt1, tt5, 0x31 ) );
	if ( i7 != B2_NULL_INDEX && ( states[i7].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i7 ), _mm256_permute2f128_ps( tt2, tt6, 0x31 ) );
	if ( i8 != B2_NULL_INDEX && ( states[i8].flags & b2_dynamicFlag ) != 0 )
		_mm256_store_ps( (float*)( states + i8 ), _mm256_permute2f128_ps( tt3, tt7, 0x31 ) );
}

#elif defined( B2_SIMD_NEON )

// This is a load and transpose
static b2BodyStateW b2GatherBodies( const b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices )
{
	_Static_assert( sizeof( b2BodyState ) == 32, "b2BodyState not 32 bytes" );
	B2_ASSERT( ( (uintptr_t)states & 0x1F ) == 0 );

	// [vx vy w flags]
	b2FloatW identityA = b2ZeroW();

	// [dpx dpy dqc dqs]

	b2FloatW identityB = b2SetW( 0.0f, 0.0f, 1.0f, 0.0f );

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;

	b2FloatW b1a = i1 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i1 ) + 0 );
	b2FloatW b1b = i1 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i1 ) + 4 );
	b2FloatW b2a = i2 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i2 ) + 0 );
	b2FloatW b2b = i2 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i2 ) + 4 );
	b2FloatW b3a = i3 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i3 ) + 0 );
	b2FloatW b3b = i3 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i3 ) + 4 );
	b2FloatW b4a = i4 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i4 ) + 0 );
	b2FloatW b4b = i4 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i4 ) + 4 );

	// [vx1 vx3 vy1 vy3]
	b2FloatW t1a = b2UnpackLoW( b1a, b3a );

	// [vx2 vx4 vy2 vy4]
	b2FloatW t2a = b2UnpackLoW( b2a, b4a );

	// [w1 w3 f1 f3]
	b2FloatW t3a = b2UnpackHiW( b1a, b3a );

	// [w2 w4 f2 f4]
	b2FloatW t4a = b2UnpackHiW( b2a, b4a );

	b2BodyStateW simdBody;
	simdBody.v.X = b2UnpackLoW( t1a, t2a );
	simdBody.v.Y = b2UnpackHiW( t1a, t2a );
	simdBody.w = b2UnpackLoW( t3a, t4a );
	simdBody.flags = b2UnpackHiW( t3a, t4a );

	b2FloatW t1b = b2UnpackLoW( b1b, b3b );
	b2FloatW t2b = b2UnpackLoW( b2b, b4b );
	b2FloatW t3b = b2UnpackHiW( b1b, b3b );
	b2FloatW t4b = b2UnpackHiW( b2b, b4b );

	simdBody.dp.X = b2UnpackLoW( t1b, t2b );
	simdBody.dp.Y = b2UnpackHiW( t1b, t2b );
	simdBody.dq.C = b2UnpackLoW( t3b, t4b );
	simdBody.dq.S = b2UnpackHiW( t3b, t4b );

	return simdBody;
}

// This writes only the velocities back to the solver bodies
// https://developer.arm.com/documentation/102107a/0100/Floating-point-4x4-matrix-transposition
static void b2ScatterBodies( b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices, const b2BodyStateW* B2_RESTRICT simdBody )
{
	_Static_assert( sizeof( b2BodyState ) == 32, "b2BodyState not 32 bytes" );
	B2_ASSERT( ( (uintptr_t)states & 0x1F ) == 0 );

	//	b2FloatW x = b2SetW(0.0f, 1.0f, 2.0f, 3.0f);
	//	b2FloatW y = b2SetW(4.0f, 5.0f, 6.0f, 7.0f);
	//	b2FloatW z = b2SetW(8.0f, 9.0f, 10.0f, 11.0f);
	//	b2FloatW w = b2SetW(12.0f, 13.0f, 14.0f, 15.0f);
	//
	//	float32x4x2_t rr1 = vtrnq_f32( x, y );
	//	float32x4x2_t rr2 = vtrnq_f32( z, w );
	//
	//	float32x4_t b1 = vcombine_f32(vget_low_f32(rr1.val[0]), vget_low_f32(rr2.val[0]));
	//	float32x4_t b2 = vcombine_f32(vget_low_f32(rr1.val[1]), vget_low_f32(rr2.val[1]));
	//	float32x4_t b3 = vcombine_f32(vget_high_f32(rr1.val[0]), vget_high_f32(rr2.val[0]));
	//	float32x4_t b4 = vcombine_f32(vget_high_f32(rr1.val[1]), vget_high_f32(rr2.val[1]));

	// transpose
	float32x4x2_t r1 = vtrnq_f32( simdBody->v.X, simdBody->v.Y );
	float32x4x2_t r2 = vtrnq_f32( simdBody->w, simdBody->flags );

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;

	// I don't use any dummy body in the body array because this will lead to multithreaded sharing and the
	// associated cache flushing.
	if ( i1 != B2_NULL_INDEX && ( states[i1].flags & b2_dynamicFlag ) != 0 )
	{
		float32x4_t body1 = vcombine_f32( vget_low_f32( r1.val[0] ), vget_low_f32( r2.val[0] ) );
		b2StoreW( (float*)( states + i1 ), body1 );
	}

	if ( i2 != B2_NULL_INDEX && ( states[i2].flags & b2_dynamicFlag ) != 0 )
	{
		float32x4_t body2 = vcombine_f32( vget_low_f32( r1.val[1] ), vget_low_f32( r2.val[1] ) );
		b2StoreW( (float*)( states + i2 ), body2 );
	}

	if ( i3 != B2_NULL_INDEX && ( states[i3].flags & b2_dynamicFlag ) != 0 )
	{
		float32x4_t body3 = vcombine_f32( vget_high_f32( r1.val[0] ), vget_high_f32( r2.val[0] ) );
		b2StoreW( (float*)( states + i3 ), body3 );
	}

	if ( i4 != B2_NULL_INDEX && ( states[i4].flags & b2_dynamicFlag ) != 0 )
	{
		float32x4_t body4 = vcombine_f32( vget_high_f32( r1.val[1] ), vget_high_f32( r2.val[1] ) );
		b2StoreW( (float*)( states + i4 ), body4 );
	}
}

#elif defined( B2_SIMD_SSE2 )

// This is a load and transpose
static b2BodyStateW b2GatherBodies( const b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices )
{
	_Static_assert( sizeof( b2BodyState ) == 32, "b2BodyState not 32 bytes" );
	B2_ASSERT( ( (uintptr_t)states & 0x1F ) == 0 );
	B2_VALIDATE( indices[0] >= 0 && indices[1] >= 0 && indices[2] >= 0 && indices[3] >= 0 );

	// [vx vy w flags]
	b2FloatW identityA = b2ZeroW();

	// [dpx dpy dqc dqs]
	b2FloatW identityB = b2SetW( 0.0f, 0.0f, 1.0f, 0.0f );

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;

	b2FloatW b1a = i1 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i1 ) + 0 );
	b2FloatW b1b = i1 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i1 ) + 4 );
	b2FloatW b2a = i2 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i2 ) + 0 );
	b2FloatW b2b = i2 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i2 ) + 4 );
	b2FloatW b3a = i3 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i3 ) + 0 );
	b2FloatW b3b = i3 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i3 ) + 4 );
	b2FloatW b4a = i4 == B2_NULL_INDEX ? identityA : b2LoadW( (float*)( states + i4 ) + 0 );
	b2FloatW b4b = i4 == B2_NULL_INDEX ? identityB : b2LoadW( (float*)( states + i4 ) + 4 );

	// [vx1 vx3 vy1 vy3]
	b2FloatW t1a = b2UnpackLoW( b1a, b3a );

	// [vx2 vx4 vy2 vy4]
	b2FloatW t2a = b2UnpackLoW( b2a, b4a );

	// [w1 w3 f1 f3]
	b2FloatW t3a = b2UnpackHiW( b1a, b3a );

	// [w2 w4 f2 f4]
	b2FloatW t4a = b2UnpackHiW( b2a, b4a );

	b2BodyStateW simdBody;
	simdBody.v.X = b2UnpackLoW( t1a, t2a );
	simdBody.v.Y = b2UnpackHiW( t1a, t2a );
	simdBody.w = b2UnpackLoW( t3a, t4a );
	simdBody.flags = b2UnpackHiW( t3a, t4a );

	b2FloatW t1b = b2UnpackLoW( b1b, b3b );
	b2FloatW t2b = b2UnpackLoW( b2b, b4b );
	b2FloatW t3b = b2UnpackHiW( b1b, b3b );
	b2FloatW t4b = b2UnpackHiW( b2b, b4b );

	simdBody.dp.X = b2UnpackLoW( t1b, t2b );
	simdBody.dp.Y = b2UnpackHiW( t1b, t2b );
	simdBody.dq.C = b2UnpackLoW( t3b, t4b );
	simdBody.dq.S = b2UnpackHiW( t3b, t4b );

	return simdBody;
}

// This writes only the velocities back to the solver bodies
static void b2ScatterBodies( b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices, const b2BodyStateW* B2_RESTRICT simdBody )
{
	_Static_assert( sizeof( b2BodyState ) == 32, "b2BodyState not 32 bytes" );
	B2_ASSERT( ( (uintptr_t)states & 0x1F ) == 0 );
	B2_VALIDATE( indices[0] >= 0 && indices[1] >= 0 && indices[2] >= 0 && indices[3] >= 0 );

	// [vx1 vy1 vx2 vy2]
	b2FloatW t1 = b2UnpackLoW( simdBody->v.X, simdBody->v.Y );
	// [vx3 vy3 vx4 vy4]
	b2FloatW t2 = b2UnpackHiW( simdBody->v.X, simdBody->v.Y );
	// [w1 f1 w2 f2]
	b2FloatW t3 = b2UnpackLoW( simdBody->w, simdBody->flags );
	// [w3 f3 w4 f4]
	b2FloatW t4 = b2UnpackHiW( simdBody->w, simdBody->flags );

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;

#if 1
	// I don't use any dummy body in the body array because this will lead to multithreaded cache coherence problems.
	if ( i1 != B2_NULL_INDEX && ( states[i1].flags & b2_dynamicFlag ) != 0 )
	{
		// [t1.x t1.y t3.x t3.y]
		b2StoreW( (float*)( states + i1 ), _mm_shuffle_ps( t1, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) ) );
	}

	if ( i2 != B2_NULL_INDEX && ( states[i2].flags & b2_dynamicFlag ) != 0 )
	{
		// [t1.z t1.w t3.z t3.w]
		b2StoreW( (float*)( states + i2 ), _mm_shuffle_ps( t1, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) ) );
	}

	if ( i3 != B2_NULL_INDEX && ( states[i3].flags & b2_dynamicFlag ) != 0 )
	{
		// [t2.x t2.y t4.x t4.y]
		b2StoreW( (float*)( states + i3 ), _mm_shuffle_ps( t2, t4, _MM_SHUFFLE( 1, 0, 1, 0 ) ) );
	}

	if ( i4 != B2_NULL_INDEX && ( states[i4].flags & b2_dynamicFlag ) != 0 )
	{
		// [t2.z t2.w t4.z t4.w]
		b2StoreW( (float*)( states + i4 ), _mm_shuffle_ps( t2, t4, _MM_SHUFFLE( 3, 2, 3, 2 ) ) );
	}

#else

	// todo_testing this is here to test the impact of unsafe writes

	if ( i1 != B2_NULL_INDEX )
	{
		// [t1.x t1.y t3.x t3.y]
		b2StoreW( (float*)( states + i1 ), _mm_shuffle_ps( t1, t3, _MM_SHUFFLE( 1, 0, 1, 0 ) ) );
	}

	if ( i2 != B2_NULL_INDEX )
	{
		// [t1.z t1.w t3.z t3.w]
		b2StoreW( (float*)( states + i2 ), _mm_shuffle_ps( t1, t3, _MM_SHUFFLE( 3, 2, 3, 2 ) ) );
	}

	if ( i3 != B2_NULL_INDEX )
	{
		// [t2.x t2.y t4.x t4.y]
		b2StoreW( (float*)( states + i3 ), _mm_shuffle_ps( t2, t4, _MM_SHUFFLE( 1, 0, 1, 0 ) ) );
	}

	if ( i4 != B2_NULL_INDEX )
	{
		// [t2.z t2.w t4.z t4.w]
		b2StoreW( (float*)( states + i4 ), _mm_shuffle_ps( t2, t4, _MM_SHUFFLE( 3, 2, 3, 2 ) ) );
	}

#endif
}

#else

// This is a load and transpose
static b2BodyStateW b2GatherBodies( const b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices )
{
	B2_VALIDATE( indices[0] >= 0 && indices[1] >= 0 && indices[2] >= 0 && indices[3] >= 0 );

	b2BodyState identity = b2_identityBodyState;

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;

	b2BodyState s1 = i1 == B2_NULL_INDEX ? identity : states[i1];
	b2BodyState s2 = i2 == B2_NULL_INDEX ? identity : states[i2];
	b2BodyState s3 = i3 == B2_NULL_INDEX ? identity : states[i3];
	b2BodyState s4 = i4 == B2_NULL_INDEX ? identity : states[i4];

	b2BodyStateW simdBody;
	simdBody.v.X = (b2FloatW){ s1.linearVelocity.x, s2.linearVelocity.x, s3.linearVelocity.x, s4.linearVelocity.x };
	simdBody.v.Y = (b2FloatW){ s1.linearVelocity.y, s2.linearVelocity.y, s3.linearVelocity.y, s4.linearVelocity.y };
	simdBody.w = (b2FloatW){ s1.angularVelocity, s2.angularVelocity, s3.angularVelocity, s4.angularVelocity };
	simdBody.flags = (b2FloatW){ (float)s1.flags, (float)s2.flags, (float)s3.flags, (float)s4.flags };
	simdBody.dp.X = (b2FloatW){ s1.deltaPosition.x, s2.deltaPosition.x, s3.deltaPosition.x, s4.deltaPosition.x };
	simdBody.dp.Y = (b2FloatW){ s1.deltaPosition.y, s2.deltaPosition.y, s3.deltaPosition.y, s4.deltaPosition.y };
	simdBody.dq.C = (b2FloatW){ s1.deltaRotation.c, s2.deltaRotation.c, s3.deltaRotation.c, s4.deltaRotation.c };
	simdBody.dq.S = (b2FloatW){ s1.deltaRotation.s, s2.deltaRotation.s, s3.deltaRotation.s, s4.deltaRotation.s };

	return simdBody;
}

// This writes only the velocities back to the solver bodies
static void b2ScatterBodies( b2BodyState* B2_RESTRICT states, int* B2_RESTRICT indices, const b2BodyStateW* B2_RESTRICT simdBody )
{
	B2_VALIDATE( indices[0] >= 0 && indices[1] >= 0 && indices[2] >= 0 && indices[3] >= 0 );

	// zero means null
	int i1 = indices[0] - 1;
	int i2 = indices[1] - 1;
	int i3 = indices[2] - 1;
	int i4 = indices[3] - 1;

	if ( i1 != B2_NULL_INDEX && ( states[i1].flags & b2_dynamicFlag ) != 0 )
	{
		b2BodyState* state = states + i1;
		state->linearVelocity.x = simdBody->v.X.x;
		state->linearVelocity.y = simdBody->v.Y.x;
		state->angularVelocity = simdBody->w.x;
	}

	if ( i2 != B2_NULL_INDEX && ( states[i2].flags & b2_dynamicFlag ) != 0 )
	{
		b2BodyState* state = states + i2;
		state->linearVelocity.x = simdBody->v.X.y;
		state->linearVelocity.y = simdBody->v.Y.y;
		state->angularVelocity = simdBody->w.y;
	}

	if ( i3 != B2_NULL_INDEX && ( states[i3].flags & b2_dynamicFlag ) != 0 )
	{
		b2BodyState* state = states + i3;
		state->linearVelocity.x = simdBody->v.X.z;
		state->linearVelocity.y = simdBody->v.Y.z;
		state->angularVelocity = simdBody->w.z;
	}

	if ( i4 != B2_NULL_INDEX && ( states[i4].flags & b2_dynamicFlag ) != 0 )
	{
		b2BodyState* state = states + i4;
		state->linearVelocity.x = simdBody->v.X.w;
		state->linearVelocity.y = simdBody->v.Y.w;
		state->angularVelocity = simdBody->w.w;
	}
}

#endif

#endif
