//=====================================================================================================
// EnhancedMadgwickAHRS.h
//=====================================================================================================
//
// Enhanced implementation of Madgwick's AHRS algorithm combining:
// - Gradient descent from original Madgwick
// - PI controller from Mahony
// - Runge-Kutta integration
//
// Date         Author          Notes
// 2025/07/04   Nanwan        Combined implementation
//
//=====================================================================================================
#ifndef EnhancedMadgwickAHRS_h
#define EnhancedMadgwickAHRS_h

//----------------------------------------------------------------------------------------------------
// Configuration
#ifndef M_PI
#define M_PI 3.14159265358979323846f
#endif
#define SAMPLE_FREQ     100.0f      // Sample frequency in Hz
#define BETA_DEF        0.1f        // Default gradient descent gain
#define KP_DEF          2.0f        // Default proportional gain
#define KI_DEF          0.005f      // Default integral gain

//----------------------------------------------------------------------------------------------------
// Variable declaration

extern volatile float beta;                     // Gradient descent gain
extern volatile float Kp, Ki;                   // PI controller gains
extern volatile float q0, q1, q2, q3;          // Quaternion elements
extern volatile float integralFBx, integralFBy, integralFBz;  // Integral error terms

//----------------------------------------------------------------------------------------------------
// Function declarations

// Main AHRS update function with magnetometer
void EnhancedMadgwickAHRSupdate(float gx, float gy, float gz, 
                                float ax, float ay, float az, 
                                float mx, float my, float mz);

// IMU update function without magnetometer
void EnhancedMadgwickAHRSupdateIMU(float gx, float gy, float gz, 
                                   float ax, float ay, float az);

// Get Euler angles from quaternion
void EnhancedMadgwickGetEulerAngles(float* roll, float* pitch, float* yaw);

// Reset quaternion to initial state
void EnhancedMadgwickReset(void);

// Set algorithm parameters
void EnhancedMadgwickSetBeta(float newBeta);
void EnhancedMadgwickSetPIGains(float newKp, float newKi);

#endif
//=====================================================================================================
// End of file
//=====================================================================================================