/*
 * Copyright (c) 2017
 *	Side Effects Software Inc.  All rights reserved.
 *
 * Redistribution and use of Houdini Development Kit samples in source and
 * binary forms, with or without modification, are permitted provided that the
 * following conditions are met:
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. The name of Side Effects Software may not be used to endorse or
 *    promote products derived from this software without specific prior
 *    written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY SIDE EFFECTS SOFTWARE `AS IS' AND ANY EXPRESS
 * OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
 * OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN
 * NO EVENT SHALL SIDE EFFECTS SOFTWARE BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA,
 * OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
 * EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *----------------------------------------------------------------------------
 */

#include "a3_rot_vec.h"
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <PRM/PRM_Include.h>
#include <SIM/SIM_PRMShared.h>
#include <SIM/SIM_DopDescription.h>
#include <SIM/SIM_ScalarField.h>
#include <SIM/SIM_VectorField.h>
#include <SIM/SIM_MatrixField.h>
#include <SIM/SIM_Object.h>
#include <GAS/GAS_SubSolver.h>

using namespace HDK_Sample;

///
/// This is the hook that Houdini grabs from the dll to link in
/// this.  As such, it merely has to implement the data factory
/// for this node.
///
void
initializeSIM(void *)
{
    IMPLEMENT_DATAFACTORY(a3RotVec);
}

/// Standard constructor, note that BaseClass was crated by the
/// DECLARE_DATAFACTORY and provides an easy way to chain through
/// the class hierarchy.
a3RotVec::a3RotVec(const SIM_DataFactory *factory)
    : BaseClass(factory)
{
}

a3RotVec::~a3RotVec()
{
}

/// Used to automatically populate the node which will represent
/// this data type.
const SIM_DopDescription *
a3RotVec::getDopDescription()
{
	static PRM_Name     theDstFieldName(GAS_NAME_FIELDDEST, "Dest Field");
	static PRM_Name     theSrcFieldName(GAS_NAME_FIELDSOURCE, "Rotor");

    static PRM_Template		 theTemplates[] = {
	PRM_Template(PRM_STRING, 1, &theDstFieldName),
	PRM_Template(PRM_STRING, 1, &theSrcFieldName),
	PRM_Template()
    };

    static SIM_DopDescription	 theDopDescription(
	    true,		
	    "a3_rot_vec",	 
	    "a3 Rotate Vec",
	    "$OS",		
	    classname(),	
	    theTemplates);	

	theDopDescription.setDefaultUniqueDataName(1);

    return &theDopDescription;
}

bool
a3RotVec::solveGasSubclass(SIM_Engine &engine,
			SIM_Object *obj,
			SIM_Time time,
			SIM_Time timestep)
{
	SIM_VectorField     *srcvector, *dstvector;
	SIM_DataArray        src, dst;

	getMatchingData(src, obj, GAS_NAME_FIELDSOURCE);
	getMatchingData(dst, obj, GAS_NAME_FIELDDEST);

	if (!dst.entries() || !src.entries())
	{
		addError(obj, SIM_MESSAGE, "Fewer fields specified.", UT_ERROR_WARNING);
		return true;
	}

	dstvector = SIM_DATA_CAST(dst(0), SIM_VectorField);
	srcvector = SIM_DATA_CAST(src(0), SIM_VectorField);

	if (dstvector && srcvector)
	{
		addFields(dstvector, srcvector);
	}

	if (dstvector)  dstvector->pubHandleModification();

	return true;

}

void
a3RotVec::addFieldsPartial(SIM_VectorField *dst, SIM_VectorField *rot, const UT_JobInfo &info)
{

	SIM_RawField *dstx = dst->getField(0);
	SIM_RawField *dsty = dst->getField(1);
	SIM_RawField *dstz = dst->getField(2);

	SIM_RawField *rotx = rot->getField(0);
	SIM_RawField *roty = rot->getField(1);
	SIM_RawField *rotz = rot->getField(2);

	UT_VoxelArrayF* dstx_rf = dstx->fieldNC();
	UT_VoxelArrayF* dsty_rf = dsty->fieldNC();
	UT_VoxelArrayF* dstz_rf = dstz->fieldNC();

	UT_VoxelArrayF* rotx_rf = rotx->fieldNC();
	UT_VoxelArrayF* roty_rf = roty->fieldNC();
	UT_VoxelArrayF* rotz_rf = rotz->fieldNC();

	UT_VoxelArrayIteratorF vit_x;

	if (!dstx->isAligned(rotx)) return;
	if (!dsty->isAligned(roty)) return;
	if (!dstz->isAligned(rotz)) return;

	int ix, iy, iz;
	int xrez, yrez, zrez;

	dstx->getVoxelRes(xrez, yrez, zrez);

	vit_x.setArray(dstx->fieldNC());
	vit_x.setPartialRange(info.job(), info.numJobs());

	UT_Vector3F vel, noise;
	UT_QuaternionF   quat;
	float            len;
	UT_Matrix3       mat;

	for (vit_x.rewind(); !vit_x.atEnd(); vit_x.advance())
	{
		if (vit_x.x() == xrez) continue;

		ix = vit_x.x();
		iy = vit_x.y();
		iz = vit_x.z();

		vel.assign(vit_x.getValue(), dsty_rf->getValue(ix, iy, iz), dstz_rf->getValue(ix, iy, iz));

		if (vel.length() > 0.1) {

			noise.assign(rotx_rf->getValue(ix, iy, iz), roty_rf->getValue(ix, iy, iz), rotz_rf->getValue(ix, iy, iz));

			if (noise.length() > 0.0)
			{
				len = noise.normalize();

				quat.updateFromAngleAxis(len, noise, 0);
				quat.getRotationMatrix(mat);

				vel *= mat;

				vit_x.setValue(vel.x());
				dsty_rf->setValue(ix, iy, iz, vel.y());
				dstz_rf->setValue(ix, iy, iz, vel.z());
			}


		}

	}
}
