// Filename : IBZebraFishkinematics3d.h

//
// Copyright (c) 2002-2019, Amneet Bhalla and Boyce Griffith
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

#ifndef included_IBZebraFishKinematics3d
#define included_IBZebraFishKinematics3d

/////////////////////////////////////// INCLUDES ///////////////////////////////////////////
#include "ibamr/ConstraintIBKinematics.h"

#include "PatchHierarchy.h"
#include "tbox/Array.h"
#include "tbox/Database.h"
#include "tbox/Pointer.h"

#include <iostream>
#include <vector>

namespace mu
{
class Parser;
}

namespace IBAMR
{
/*!
 * \brief IBZebraFishKinematics3d Class.
 *
 * IBZebraFishKinematics3d is a concrete class which calculates the deformation velocity and updated shape
 * for a 3d Zebrafish.
 */

class IBZebraFishKinematics3d : public IBAMR::ConstraintIBKinematics
{
public:
    /*!
     * \brief Constructor.
     */
    IBZebraFishKinematics3d(const std::string& object_name,
                      SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> input_db,
                      IBTK::LDataManager* l_data_manager,
                      SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy,
                      bool register_for_restart = true);

    /*!
     * \brief Destructor.
     */
    virtual ~IBZebraFishKinematics3d();

    /*!
     * \brief Set kinematics velocity for the ed eel at specified time.
     * \see IBAMR::ConstraintIBKinematics::setKinematicsVelocity
     */
    virtual void setKinematicsVelocity(const double time,
                                       const std::vector<double>& incremented_angle_from_reference_axis,
                                       const std::vector<double>& center_of_mass,
                                       const std::vector<double>& tagged_pt_position);

    /*!
     * \brief Get the kinematics velocity on the specified level.
     * \see IBAMR::ConstraintIBKinematics::getKinematicsVelocity
     */
    virtual const std::vector<std::vector<double> >& getKinematicsVelocity(const int level) const;

    /*!
     * \brief Set the shape of eel at specified time. The shape should have
     * its center of mass at the origin, with appropriate rigid body rotation applied
     * to it.
     * \see IBAMR::ConstraintIBKinematics::setShape
     */
    virtual void setShape(const double time, const std::vector<double>& incremented_angle_from_reference_axis);

    /*!
     * \brief Get the shape of eel on the specified level. The shape should have
     * its center of mass at the origin, with appropriate rigid body rotation applied
     * to it.
     * \see IBAMR::ConstraintIBKinematics::getShape
     */
    virtual const std::vector<std::vector<double> >& getShape(const int level) const;

    /*!
     * \brief Write state necessary for restarted runs.
     */
    virtual void putToDatabase(SAMRAI::tbox::Pointer<SAMRAI::tbox::Database> db);

private:
    /*!
     * \brief The default constructor is not implemented and should not be used.
     */
    IBZebraFishKinematics3d();

    /*!
     * \brief The copy constructor is not implemented and should not be used.
     */
    IBZebraFishKinematics3d(const IBZebraFishKinematics3d& from);

    /*!
     * \brief The assignment operator is not implemented and should not be used.
     */
    IBZebraFishKinematics3d& operator=(const IBZebraFishKinematics3d& that);

    /*!
     * \brief Calculate the coefficients of the cubic polynomial
     * from the curvature coefficients.
     */
    void getInterpCoefs();

    /*!
     * \brief Set data from restart.
     */
    void getFromRestart();

    /*!
     * \brief Set eel body shape related data.
     */
    void setImmersedBodyLayout(SAMRAI::tbox::Pointer<SAMRAI::hier::PatchHierarchy<NDIM> > patch_hierarchy);

    /*!
     * \brief Set eel kinematics velocity.
     */
    void setZebraFishSpecificVelocity(const double time,
                                const std::vector<double>& incremented_angle_from_reference_axis,
                                const std::vector<double>& center_of_mass,
                                const std::vector<double>& tagged_pt_position);

    /*!
     * Current time (t) and new time (t+dt).
     */
    double d_current_time, d_new_time;

    /*!
     * Deformational velocity and shape vectors.
     */
    std::vector<std::vector<double> > d_kinematics_vel;
    std::vector<std::vector<double> > d_shape;

    /*!
     * Save COM, tagged point position and incremented angle from reference axis for restarted runs.
     */
    std::vector<double> d_center_of_mass, d_incremented_angle_from_reference_axis, d_tagged_pt_position;

    /*!
     * Time period of traveling wave.
     */
    double d_time_period;

    /*!
     * Initial angle of the body axis from the horizontal.
     */
    double d_initAngle_bodyAxis_x;

    /*!
     * No. of points on eel.
     */
    int d_BodyNs;



    /*!
     * maximum timestep size
     */
    double d_max_time_step;

    /*!
     * Number of points composing the backbone of the 2D recording
     */
    int d_num_frames;

    int d_num_pts_backbone;

    int d_iteration_num;

    int d_num_points_model;

    int num_pts;



    /*!
     * Inputting 2D backbone coordinates
     */
    int row_counter;

    int column_counter;

    int time_ticker;

    double datum; 

    int datumm;

    double datummm;

    double datummmm;

    double datummmmm;


    std::vector<std::vector<double> > d_slice_sizes;

    std::vector<std::vector<double> > d_dummy_model;

    std::vector<std::vector<double> > heading_angle;


    /*!
     * Storing 2D backbone coordinates
     */
    std::vector<std::vector<std::vector<double> > > d_backbone_spline;




    /*!
     * Mapping 2D to 3D
     */
    std::vector<double> d_xval;
    std::vector<double> d_yval;
    int initial_tframe, x_index, y_index, z_index, init_slice, cumu_col;
    

}; // IBEELKinematics3d

} // namespace IBAMR

#endif //#ifndef included_IBEELKinematics3d