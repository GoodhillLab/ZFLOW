// Filename: IBZebraFishKinematics3d.cpp
//
// Copyright (c) 2002-2014, Amneet Bhalla and Boyce Griffith
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//    * Redistributions of source code must retain the above copyright notice,
//      this list of conditions and the following disclaimer.
//
//    * Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//    * Neither the name of The University of North Carolina nor the names of
//      its contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
// LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
// CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
// SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
// INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
// CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
// ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

//////////////////////////////////// INCLUDES ////////////////////////////////////////////

#include "ibamr/namespaces.h"
#include "CartesianPatchGeometry.h"
#include "IBZebraFishKinematics3d.h"
#include "PatchLevel.h"
#include "tbox/MathUtilities.h"
#include "tbox/SAMRAI_MPI.h"
#include "tbox/Utilities.h"
#include "muParser.h"
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <algorithm>

namespace IBAMR
{
IBZebraFishKinematics3d::IBZebraFishKinematics3d(const std::string& object_name,
                                     Pointer<Database> input_db,
                                     LDataManager* l_data_manager,
                                     Pointer<PatchHierarchy<NDIM> > patch_hierarchy,
                                     bool register_for_restart)
        : ConstraintIBKinematics(object_name, input_db, l_data_manager, register_for_restart),
          d_current_time(0.0),
          d_kinematics_vel(NDIM),
          d_shape(NDIM),
          d_center_of_mass(3), //There are n=3 movement axis (x,y,z planes)
          d_incremented_angle_from_reference_axis(3),
          d_tagged_pt_position(3),
          d_num_frames(),
          d_num_pts_backbone(),
          d_iteration_num(0),
          row_counter(),
          column_counter(),
          time_ticker(),
          d_slice_sizes(NDIM,std::vector<double>(125)), //There are n=125 points describing the 2D larval zebrafish movement
          d_dummy_model(NDIM),
          heading_angle(125),
          d_backbone_spline(125,std::vector<std::vector<double> > (2,std::vector<double>()))
{
    // Input solver details
    d_initAngle_bodyAxis_x = input_db->getDoubleWithDefault("initial_angle_horizontal", 0.0);
    d_max_time_step = input_db->getDouble("time_step_size");
    d_num_frames = input_db->getInteger("number_frames");
    d_num_pts_backbone = input_db->getInteger("num_backbone_points");
    d_num_points_model = input_db->getInteger("num_points_model");

    // Input 2D/3D data
    double datum, datumm, datummm, datummmmm;
    // READ THE TRACKING DATA
    // Resize the innermost vectors
    for (auto& outerVec : d_backbone_spline) {
        for (auto& middleVec : outerVec) {
            middleVec.resize(d_num_frames); 
        }
    }
    std::ifstream input_file("zebrafish_2D_coords.txt"); // create and open
    row_counter = 0;
    column_counter = 0;
    time_ticker = 0;  
    // Read and store input values
    while (input_file >> datum) {
        d_backbone_spline[row_counter][column_counter][time_ticker]  = datum;
        if (column_counter == 1) {
            column_counter = -1;
            row_counter++;
            //if end of column, drop down and input the next row
            if (row_counter == d_num_pts_backbone) {
                time_ticker++;
                row_counter = 0;
                //If final coordinate for a time-period, move onto next time period. Reset all other counters
            }
        }
        //Update column index
        column_counter++;
    }
    input_file.close();


    // READ THE SLICE INFORMATION
    std::ifstream input_fileee("sizes.txt"); // create and open
    row_counter = 0;
    column_counter = 0;
    while (input_fileee >> datumm) {
        d_slice_sizes[column_counter][row_counter]  = datumm;
        if (column_counter == 2) {
            column_counter = -1;
            row_counter++; 
        }
        column_counter++;     
    }
    input_fileee.close();

    // READ THE 3D MODEL
    for (int d = 0; d < NDIM; ++d)
    {
        d_dummy_model[d].resize(d_num_points_model);
    }
    std::ifstream input_filee("stacked.txt"); // create and open
    row_counter = 0;
    column_counter = 0;
    // Read and store the stacked file
    while (input_filee >> datummm) {
        d_dummy_model[column_counter][row_counter] = datummm;
        if (column_counter == 2) {
            column_counter = -1;
            row_counter++;  
        }
        //Update column index
        column_counter++;
    }
    input_filee.close();


    // READ THE HEADING ANGLES
    int storing = 124*d_num_frames;
    for (int d = 0; d < 125; ++d)
    {
        heading_angle[d].resize(storing);
    }
    // Read the heading angles
    std::ifstream input_angle("heading_angle.txt"); // create and open
    row_counter = 0;
    column_counter = 0;
    while (input_angle >> datummmmm) {
        heading_angle[column_counter][row_counter] = datummmmm;
        if (column_counter == 123) {
            column_counter = -1;
            row_counter++;  
        }
        //Update column index
        column_counter++;
    }
    input_angle.close();
    

    // set how the immersed body is layout in reference frame.
    setImmersedBodyLayout(patch_hierarchy);
    bool from_restart = RestartManager::getManager()->isFromRestart();
    if (from_restart) getFromRestart();
    return;   
}


IBZebraFishKinematics3d::~IBZebraFishKinematics3d()
{
    // intentionally left blank
    return;

} // ~IBZebraFishKinematics3d


void
IBZebraFishKinematics3d::putToDatabase(Pointer<Database> db)
{
    db->putDouble("d_current_time", d_current_time);
    db->putDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->putDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->putDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // putToDatabase


void
IBZebraFishKinematics3d::getFromRestart()
{
    Pointer<Database> restart_db = RestartManager::getManager()->getRootDatabase();
    Pointer<Database> db;
    if (restart_db->isDatabase(d_object_name))
    {
        db = restart_db->getDatabase(d_object_name);
    }
    else
    {
        TBOX_ERROR(d_object_name << ":  Restart database corresponding to " << d_object_name
                                 << " not found in restart file." << std::endl);
    }

    d_current_time = db->getDouble("d_current_time");
    db->getDoubleArray("d_center_of_mass", &d_center_of_mass[0], 3);
    db->getDoubleArray("d_incremented_angle_from_reference_axis", &d_incremented_angle_from_reference_axis[0], 3);
    db->getDoubleArray("d_tagged_pt_position", &d_tagged_pt_position[0], 3);

    return;

} // getFromRestart


void
IBZebraFishKinematics3d::setImmersedBodyLayout(Pointer<PatchHierarchy<NDIM> > patch_hierarchy)
{
    const StructureParameters& struct_param = getStructureParameters();
    const int coarsest_ln = struct_param.getCoarsestLevelNumber();
    const int finest_ln = struct_param.getFinestLevelNumber();
    TBOX_ASSERT(coarsest_ln == finest_ln);
    const std::vector<std::pair<int, int> >& idx_range = struct_param.getLagIdxRange();
    const int total_lag_pts = idx_range[0].second - idx_range[0].first;
    for (int d = 0; d < NDIM; ++d)
    {
        d_kinematics_vel[d].resize(total_lag_pts);
        d_shape[d].resize(total_lag_pts);
    }
    // Get Background mesh related data.
    Pointer<PatchLevel<NDIM> > level = patch_hierarchy->getPatchLevel(coarsest_ln);
    PatchLevel<NDIM>::Iterator p(level);
    Pointer<Patch<NDIM> > patch = level->getPatch(p());
    Pointer<CartesianPatchGeometry<NDIM> > pgeom = patch->getPatchGeometry();
    //const double* const dx = pgeom->getDx(); This is is useful line, if you want to double check the cell widths for innermost layer
    //Set-up
    d_BodyNs = 125;
    initial_tframe = 0;
    x_index = 0;
    y_index = 1;
    z_index = 2;
    init_slice = 0; //Index of the first slice
    cumu_col = 0;  //Column in sizes.txt that contains the individual slice sizes
    //Store backbone position at t=0
    d_xval.resize(d_BodyNs);
    d_yval.resize(d_BodyNs);
    std::vector<double> init_shiftx(d_BodyNs);
    std::vector<double> init_shifty(d_BodyNs);
    for (int i=1; i<=d_BodyNs; ++i)
    {
        d_xval[i-1] = d_backbone_spline[i-1][x_index][initial_tframe]; 
        d_yval[i-1] = d_backbone_spline[i-1][y_index][initial_tframe];
    } 
    return;
} // setImmersedBodyLayout


void
IBZebraFishKinematics3d::setZebraFishSpecificVelocity(const double time,
                                          const std::vector<double>& incremented_angle_from_reference_axis,
                                          const std::vector<double>& /*center_of_mass*/,
                                          const std::vector<double>& /*tagged_pt_position*/)
{
    const double angleFromHorizontal =  d_initAngle_bodyAxis_x + incremented_angle_from_reference_axis[2];
    std::vector<double> vec_vel(NDIM, 0.0);
    int lag_idx = 0;
    double dydt, dxdt;    
    int current_frame = d_iteration_num;  
    for (int i = 1; i <= d_BodyNs; ++i)
    {
        //Calculate deformational velocity
        if (current_frame == 0 || current_frame == d_num_frames-1)
        {
            dxdt = 0;
            dydt = 0;
        } else {
            dxdt = (d_backbone_spline[i-1][0][current_frame + 1] - d_backbone_spline[i-1][0][current_frame - 1])/(2*d_max_time_step);
            dydt = (d_backbone_spline[i-1][1][current_frame + 1] - d_backbone_spline[i-1][1][current_frame - 1])/(2*d_max_time_step);
        }
        //Convert defVel from relative reference frame to laboratory reference frame
        vec_vel[0] = dxdt * (std::cos(angleFromHorizontal)) + dydt * (-std::sin(angleFromHorizontal));
        vec_vel[1] = dydt * (std::cos(angleFromHorizontal)) + dxdt * (std::sin(angleFromHorizontal));
        //Apply defVel to all points in given 3D slice
        int pts_this_xsection_new = d_slice_sizes[0][i-1];
        const int lowerlimit = lag_idx;
        const int upperlimit = lag_idx + pts_this_xsection_new;
        for (int d = 0; d < NDIM; ++d)
            for (int j = lowerlimit; j < upperlimit; ++j) d_kinematics_vel[d][j] = vec_vel[d];
        lag_idx = upperlimit;
    }
    return;
} // setEelSpecificvelocity


void
IBZebraFishKinematics3d::setKinematicsVelocity(const double new_time,
                                         const std::vector<double>& incremented_angle_from_reference_axis,
                                         const std::vector<double>& center_of_mass,
                                         const std::vector<double>& tagged_pt_position)
{
    d_new_time = new_time;
    d_incremented_angle_from_reference_axis = incremented_angle_from_reference_axis;
    d_center_of_mass = center_of_mass;
    d_tagged_pt_position = tagged_pt_position;
    setZebraFishSpecificVelocity(d_new_time, d_incremented_angle_from_reference_axis, d_center_of_mass, d_tagged_pt_position);
    return;

} // setKinematicsVelocity


const std::vector<std::vector<double> >&
IBZebraFishKinematics3d::getKinematicsVelocity(const int /*level*/) const
{
    return d_kinematics_vel;

} // getNewKinematicsVelocity


void
IBZebraFishKinematics3d::setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis)
{
    const StructureParameters& struct_param = getStructureParameters();
    const std::string position_update_method = struct_param.getPositionUpdateMethod();
    if (position_update_method == "CONSTRAINT_VELOCITY")
    {
        return;
    }
    else
    {        
        //Set-up
        int current_frame = d_iteration_num;
        TBOX_ASSERT(d_new_time == time);
        double xbase, ybase;
        std::vector<double> shape_new(NDIM);
        int lag_idx = -1;
        //Reconstruct the 3D fish around the 2D backbone
        int start = 0;
        int finish = start + d_slice_sizes[cumu_col][init_slice]-1;
        for (int i = 1; i <= d_BodyNs; ++i)
        {
            //
            // PART ONE: MAPPING BETWEEN 2D TAILTRACKING AND 3D MODEL::
            //
            //Grab position of the 2D recording
            xbase = d_backbone_spline[i-1][x_index][current_frame];
            ybase = d_backbone_spline[i-1][y_index][current_frame];	
            // Calculate the total shift of dummy eel3d.vertex
            double shift_x = xbase - d_xval[i-1];    
            double shift_y = ybase - d_yval[i-1];
            // Temporary matrix to hold slice
            std::vector<std::vector<double> > spin_slice(NDIM);
            for (int d = 0; d < NDIM; ++d)
            {
                spin_slice[d].resize(finish-start+1);
            }
            // Map slice onto 2D movement
            int cc = 0;
            for (int j = start; j <= finish; ++j)
            {
                spin_slice[x_index][cc] = d_dummy_model[x_index][j] + shift_x;
                spin_slice[y_index][cc] = d_dummy_model[y_index][j] + shift_y;
                spin_slice[z_index][cc] = d_dummy_model[z_index][j];
                ++cc;
            }            
            //
            // PART TWO: LOCAL SLICE ROTATION:: Here we are ensuring the shape of the fish is stable when turning (ie. that the mesh points don't mix unrealistically)
            //
            //Grab the rotation angle:
            double theta = heading_angle[i-1][current_frame]; //For debugging, or to see what happens when the below code section isn't used, set theta=0
            //Calculate centre of mass of slice:
            std::vector<double> V_Centre(NDIM);
            for (int kk = 0; kk < NDIM; ++kk) 
            {
                double hold = 0;
                for (int kkk = 0; kkk < (finish-start+1); ++kkk)
                {
                    hold = hold + spin_slice[kk][kkk];
                }
                V_Centre[kk] = hold/(finish-start+1);
            }
            //Subtract COM from the points on the slice (so that we can rotate around the origin)
            std::vector<std::vector<double> > Vc(NDIM);
            for (int d = 0; d < NDIM; ++d) {
                Vc[d].resize(finish-start+1);
            }
            cc = 0;
            for (int j = start; j <= finish; ++j)
            {
                Vc[x_index][cc] = spin_slice[x_index][cc] - V_Centre[x_index];
                Vc[y_index][cc] = spin_slice[y_index][cc] - V_Centre[y_index];
                Vc[z_index][cc] = spin_slice[z_index][cc] - V_Centre[z_index];                
                ++cc;
            }
            //Convert angle to radians
            double a_rad = ((theta*M_PI)/180);
            //Create angles matrix
            std::vector<double> angles_mat{0, 0, a_rad};
            //create rotation matrix
            std::vector<std::vector<double> > rot_mat{
                {cos(angles_mat[2])*cos(angles_mat[1]),     (cos(angles_mat[2])*sin(angles_mat[1])*sin(angles_mat[0])) - (sin(angles_mat[2])*cos(angles_mat[0])),   (cos(angles_mat[2])*sin(angles_mat[1])*cos(angles_mat[0])) + (sin(angles_mat[2])*sin(angles_mat[0])) },
                {sin(angles_mat[2])*cos(angles_mat[1]),     (sin(angles_mat[2])*sin(angles_mat[1])*sin(angles_mat[0])) + (cos(angles_mat[2])*cos(angles_mat[0])),   (sin(angles_mat[2])*sin(angles_mat[1])*cos(angles_mat[0])) - (cos(angles_mat[2])*sin(angles_mat[0])) },
                {-sin(angles_mat[1])                  ,     cos(angles_mat[1])*sin(angles_mat[0])                                                               ,   cos(angles_mat[1])*cos(angles_mat[0]) }
            };
            //Multiply each point in slice by the rotation matrix
            int r1 = 3;
            int c1 = 3;
            int c2 = finish-start+1;
            std::vector<std::vector<double> > mult(NDIM);
            for (int d = 0; d < NDIM; ++d)
            {
                mult[d].resize(c2);
            }
            // Multiplying matrix rot_mat and Vc and storing in array mult.
            for(int jj = 0; jj < r1; ++jj)
                for(int jjj = 0; jjj < c2; ++jjj)
                    for(int k = 0; k < c1; ++k)
                    {
                        mult[jj][jjj] += rot_mat[jj][k] * Vc[k][jjj];
                    }

            //Then, shift the rotated slice back to the original Centre of Mass 
            std::vector<std::vector<double> > V_r(NDIM);
            for (int d = 0; d < NDIM; ++d) {
                V_r[d].resize(finish-start+1);
            }
            cc = 0;
            for (int j = start; j <= finish; ++j)
            {
                V_r[x_index][cc] = mult[x_index][cc] + V_Centre[x_index];
                V_r[y_index][cc] = mult[y_index][cc] + V_Centre[y_index];
                V_r[z_index][cc] = mult[z_index][cc] + V_Centre[z_index];
                ++cc;
            }
            //Finally, add slice to d_shape
            cc = 0;
            for (int j = start; j <= finish; ++j)
            {
                d_shape[x_index][++lag_idx] = V_r[x_index][cc];
                d_shape[y_index][lag_idx] = V_r[y_index][cc];
                d_shape[z_index][lag_idx] = V_r[z_index][cc];
                ++cc;
            }
            //Update slice
            start = finish+1;
            if (i<125){
                finish=start+d_slice_sizes[cumu_col][i]-1; //i.e grab the next slice's NumPoints
            }
	    }
        //
        // PART THREE: NORMALIZATION::
        //
        // Calculate the C.O.M of the reconstructed 3D fish
        std::vector<double> center_of_mass(NDIM, 0.0);
        const int total_lag_pts = d_shape[0].size();
        for (int d = 0; d < NDIM; ++d)
        {
            for (std::vector<double>::const_iterator citr = d_shape[d].begin(); citr != d_shape[d].end(); ++citr)
            {
                center_of_mass[d] += *citr;
            }
        }
        for (int d = 0; d < NDIM; ++d) center_of_mass[d] /= total_lag_pts;
        // Shift the c.m to the origin to apply the rotation
        for (int d = 0; d < NDIM; ++d)
        {
            for (std::vector<double>::iterator itr = d_shape[d].begin(); itr != d_shape[d].end(); ++itr)
            {
                *itr -= center_of_mass[d];
            }
        }
        // Apply rotation
        const double angleFromHorizontal = d_initAngle_bodyAxis_x + incremented_angle_from_reference_axis[2];
        for (int i = 0; i < total_lag_pts; ++i)
        {
            const double x_rotated = d_shape[0][i] * cos(angleFromHorizontal) - d_shape[1][i] * sin(angleFromHorizontal);
            const double y_rotated = d_shape[0][i] * sin(angleFromHorizontal) + d_shape[1][i] * cos(angleFromHorizontal);
            d_shape[0][i] = x_rotated;
            d_shape[1][i] = y_rotated;
        }
        // Update time
        if (time > 0) 
        {
            ++d_iteration_num;
        }
        return;
    }
    d_current_time = d_new_time;
} // setShape


const std::vector<std::vector<double> >&
IBZebraFishKinematics3d::getShape(const int /*level*/) const
{
    return d_shape;
} // getShape

} // namespace IBAMR
