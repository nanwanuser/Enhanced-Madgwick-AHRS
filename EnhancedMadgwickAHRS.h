//=====================================================================================================
// EnhancedMadgwickAHRS.h
//=====================================================================================================
//
// Enhanced implementation of Madgwick's IMU and AHRS algorithms.
// Based on the original Madgwick algorithm with 4th-order Runge-Kutta integration.
// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
//
// Date			Author          Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
// 04/07/2025	Nanwan          Added 4th-order Runge-Kutta integration
//
//=====================================================================================================
#ifndef EnhancedMadgwickAHRS_h
#define EnhancedMadgwickAHRS_h

//----------------------------------------------------------------------------------------------------
// Variable declaration

extern volatile float beta;				// algorithm gain
extern volatile float q0, q1, q2, q3;	// quaternion of sensor frame relative to auxiliary frame

//---------------------------------------------------------------------------------------------------
// Function declarations

void EnhancedMadgwickAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);
void EnhancedMadgwickAHRSupdateIMU(float gx, float gy, float gz, float ax, float ay, float az);

#endif
//=====================================================================================================
// End of file
//=====================================================================================================