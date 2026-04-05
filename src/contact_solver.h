// SPDX-FileCopyrightText: 2023 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "solver.h"

typedef struct b2BodySim b2BodySim;
typedef struct b2ContactSim b2ContactSim;

typedef struct b2ContactConstraintPoint
{
	b2Vec2 anchorA, anchorB;
	float baseSeparation;
	float relativeVelocity;
	float normalImpulse;
	float tangentImpulse;
	float totalNormalImpulse;
	float normalMass;
	float tangentMass;
} b2ContactConstraintPoint;

typedef struct b2ContactConstraint
{
	// base-1, 0 for null
	int indexA;
	int indexB;
	b2ContactConstraintPoint points[2];
	b2Vec2 normal;
	float invMassA, invMassB;
	float invIA, invIB;
	float friction;
	float restitution;
	float tangentSpeed;
	float rollingResistance;
	float rollingMass;
	float rollingImpulse;
	b2Softness softness;
	int pointCount;
} b2ContactConstraint;

int b2GetContactConstraintSIMDByteCount( void );

// Scalar constraint functions for cluster solver (operate on arbitrary constraint arrays)
void b2PrepareContactConstraints( b2ContactSim** contacts, b2ContactConstraint* constraints, int count,
								  b2StepContext* context, b2BodySim* bodySims );
void b2WarmStartContactConstraints( b2ContactConstraint* constraints, int count, b2BodyState* states );
void b2SolveContactConstraints( b2ContactConstraint* constraints, int count, b2BodyState* states, float inv_h,
								float contactSpeed, bool useBias );
void b2ApplyContactRestitution( b2ContactConstraint* constraints, int count, b2BodyState* states, float threshold );
void b2StoreContactImpulses( b2ContactSim** contacts, b2ContactConstraint* constraints, int count );
