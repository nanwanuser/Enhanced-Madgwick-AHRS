//=====================================================================================================
// EnhancedMadgwickAHRS.c
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

//---------------------------------------------------------------------------------------------------
// Header files

#include "EnhancedMadgwickAHRS.h"
#include <math.h>

//---------------------------------------------------------------------------------------------------
// Definitions

#define dt (1.0f / SAMPLE_FREQ)    // Sample period

//---------------------------------------------------------------------------------------------------
// Variable definitions

volatile float beta = BETA_DEF;                                    // Gradient descent gain
volatile float Kp = KP_DEF;                                       // Proportional gain
volatile float Ki = KI_DEF;                                       // Integral gain
volatile float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;      // Quaternion
volatile float integralFBx = 0.0f, integralFBy = 0.0f, integralFBz = 0.0f;  // Integral error

// Runge-Kutta intermediate variables
static float k1[4], k2[4], k3[4], k4[4];

//---------------------------------------------------------------------------------------------------
// Function declarations

float invSqrt(float x);
void quaternionDerivative(float q0, float q1, float q2, float q3, 
                         float wx, float wy, float wz, float* dq);

//====================================================================================================
// Functions

//---------------------------------------------------------------------------------------------------
// Enhanced AHRS algorithm update

void EnhancedMadgwickAHRSupdate(float gx, float gy, float gz, 
                               float ax, float ay, float az, 
                               float mx, float my, float mz) {
    float recipNorm;
    float s0, s1, s2, s3;
    float qDot1, qDot2, qDot3, qDot4;
    float hx, hy;
    float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz;
    float _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3;
    float q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;
    float ex, ey, ez;
    float qa, qb, qc, qd;
    
    // Use IMU algorithm if magnetometer measurement invalid
    if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
        EnhancedMadgwickAHRSupdateIMU(gx, gy, gz, ax, ay, az);
        return;
    }
    
    // Compute feedback only if accelerometer measurement valid
    if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {
        
        // Normalize accelerometer measurement
        recipNorm = invSqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;
        
        // Normalize magnetometer measurement
        recipNorm = invSqrt(mx * mx + my * my + mz * mz);
        mx *= recipNorm;
        my *= recipNorm;
        mz *= recipNorm;
        
        // Auxiliary variables to avoid repeated arithmetic
        _2q0mx = 2.0f * q0 * mx;
        _2q0my = 2.0f * q0 * my;
        _2q0mz = 2.0f * q0 * mz;
        _2q1mx = 2.0f * q1 * mx;
        _2q0 = 2.0f * q0;
        _2q1 = 2.0f * q1;
        _2q2 = 2.0f * q2;
        _2q3 = 2.0f * q3;
        _2q0q2 = 2.0f * q0 * q2;
        _2q2q3 = 2.0f * q2 * q3;
        q0q0 = q0 * q0;
        q0q1 = q0 * q1;
        q0q2 = q0 * q2;
        q0q3 = q0 * q3;
        q1q1 = q1 * q1;
        q1q2 = q1 * q2;
        q1q3 = q1 * q3;
        q2q2 = q2 * q2;
        q2q3 = q2 * q3;
        q3q3 = q3 * q3;
        
        // Reference direction of Earth's magnetic field
        hx = mx * q0q0 - _2q0my * q3 + _2q0mz * q2 + mx * q1q1 + _2q1 * my * q2 + _2q1 * mz * q3 - mx * q2q2 - mx * q3q3;
        hy = _2q0mx * q3 + my * q0q0 - _2q0mz * q1 + _2q1mx * q2 - my * q1q1 + my * q2q2 + _2q2 * mz * q3 - my * q3q3;
        _2bx = sqrtf(hx * hx + hy * hy);
        _2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3;
        _4bx = 2.0f * _2bx;
        _4bz = 2.0f * _2bz;
        
        // Gradient descent algorithm corrective step
        s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q1 * (1.0f - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q2 * (1.0f - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
        recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;
        
        // Compute error for PI controller (cross product between estimated and measured direction)
        float vx = 2.0f * (q1q3 - q0q2);
        float vy = 2.0f * (q0q1 + q2q3);
        float vz = q0q0 - q1q1 - q2q2 + q3q3;
        
        float wx = 2.0f * _2bx * (0.5f - q2q2 - q3q3) + 2.0f * _2bz * (q1q3 - q0q2);
        float wy = 2.0f * _2bx * (q1q2 - q0q3) + 2.0f * _2bz * (q0q1 + q2q3);
        float wz = 2.0f * _2bx * (q0q2 + q1q3) + 2.0f * _2bz * (0.5f - q1q1 - q2q2);
        
        ex = (ay * vz - az * vy) + (my * wz - mz * wy);
        ey = (az * vx - ax * vz) + (mz * wx - mx * wz);
        ez = (ax * vy - ay * vx) + (mx * wy - my * wx);
        
        // Apply integral feedback
        if(Ki > 0.0f) {
            integralFBx += Ki * ex * dt;
            integralFBy += Ki * ey * dt;
            integralFBz += Ki * ez * dt;
            gx += integralFBx;
            gy += integralFBy;
            gz += integralFBz;
        } else {
            integralFBx = 0.0f;
            integralFBy = 0.0f;
            integralFBz = 0.0f;
        }
        
        // Apply proportional feedback and gradient descent  // 应该考虑是否需要权重分配：待定……
        gx += Kp * ex - beta * s1;
        gy += Kp * ey - beta * s2;
        gz += Kp * ez - beta * s3;
        
        // Store current quaternion
        qa = q0;
        qb = q1;
        qc = q2;
        qd = q3;
        
        // Runge-Kutta 4th order integration
        // k1
        quaternionDerivative(qa, qb, qc, qd, gx, gy, gz, k1);
        
        // k2
        quaternionDerivative(qa + 0.5f * dt * k1[0], 
                           qb + 0.5f * dt * k1[1],
                           qc + 0.5f * dt * k1[2],
                           qd + 0.5f * dt * k1[3],
                           gx, gy, gz, k2);
        
        // k3
        quaternionDerivative(qa + 0.5f * dt * k2[0],
                           qb + 0.5f * dt * k2[1],
                           qc + 0.5f * dt * k2[2],
                           qd + 0.5f * dt * k2[3],
                           gx, gy, gz, k3);
        
        // k4
        quaternionDerivative(qa + dt * k3[0],
                           qb + dt * k3[1],
                           qc + dt * k3[2],
                           qd + dt * k3[3],
                           gx, gy, gz, k4);
        
        // Update quaternion
        q0 += dt / 6.0f * (k1[0] + 2.0f * k2[0] + 2.0f * k3[0] + k4[0]);
        q1 += dt / 6.0f * (k1[1] + 2.0f * k2[1] + 2.0f * k3[1] + k4[1]);
        q2 += dt / 6.0f * (k1[2] + 2.0f * k2[2] + 2.0f * k3[2] + k4[2]);
        q3 += dt / 6.0f * (k1[3] + 2.0f * k2[3] + 2.0f * k3[3] + k4[3]);
        
        // Normalize quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        q0 *= recipNorm;
        q1 *= recipNorm;
        q2 *= recipNorm;
        q3 *= recipNorm;
    }
}

//---------------------------------------------------------------------------------------------------
// IMU algorithm update (no magnetometer)

void EnhancedMadgwickAHRSupdateIMU(float gx, float gy, float gz, 
                                  float ax, float ay, float az) {
    float recipNorm;
    float s0, s1, s2, s3;
    float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2, _8q1, _8q2;
    float q0q0, q1q1, q2q2, q3q3;
    float ex, ey, ez;
    float qa, qb, qc, qd;
    
    // Compute feedback only if accelerometer measurement valid
    if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {
        
        // Normalize accelerometer measurement
        recipNorm = invSqrt(ax * ax + ay * ay + az * az);
        ax *= recipNorm;
        ay *= recipNorm;
        az *= recipNorm;
        
        // Auxiliary variables
        _2q0 = 2.0f * q0;
        _2q1 = 2.0f * q1;
        _2q2 = 2.0f * q2;
        _2q3 = 2.0f * q3;
        _4q0 = 4.0f * q0;
        _4q1 = 4.0f * q1;
        _4q2 = 4.0f * q2;
        _8q1 = 8.0f * q1;
        _8q2 = 8.0f * q2;
        q0q0 = q0 * q0;
        q1q1 = q1 * q1;
        q2q2 = q2 * q2;
        q3q3 = q3 * q3;
        
        // Gradient descent algorithm corrective step
        s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
        s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
        s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
        s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
        recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3);
        s0 *= recipNorm;
        s1 *= recipNorm;
        s2 *= recipNorm;
        s3 *= recipNorm;
        
        // Compute error for PI controller
        float vx = 2.0f * (q1 * q3 - q0 * q2);
        float vy = 2.0f * (q0 * q1 + q2 * q3);
        float vz = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
        
        ex = (ay * vz - az * vy);
        ey = (az * vx - ax * vz);
        ez = (ax * vy - ay * vx);
        
        // Apply integral feedback
        if(Ki > 0.0f) {
            integralFBx += Ki * ex * dt;
            integralFBy += Ki * ey * dt;
            integralFBz += Ki * ez * dt;
            gx += integralFBx;
            gy += integralFBy;
            gz += integralFBz;
        } else {
            integralFBx = 0.0f;
            integralFBy = 0.0f;
            integralFBz = 0.0f;
        }
        
        // Apply proportional feedback and gradient descent
        gx += Kp * ex - beta * s1;
        gy += Kp * ey - beta * s2;
        gz += Kp * ez - beta * s3;
        
        // Store current quaternion
        qa = q0;
        qb = q1;
        qc = q2;
        qd = q3;
        
        // Runge-Kutta 4th order integration
        // k1
        quaternionDerivative(qa, qb, qc, qd, gx, gy, gz, k1);
        
        // k2
        quaternionDerivative(qa + 0.5f * dt * k1[0], 
                           qb + 0.5f * dt * k1[1],
                           qc + 0.5f * dt * k1[2],
                           qd + 0.5f * dt * k1[3],
                           gx, gy, gz, k2);
        
        // k3
        quaternionDerivative(qa + 0.5f * dt * k2[0],
                           qb + 0.5f * dt * k2[1],
                           qc + 0.5f * dt * k2[2],
                           qd + 0.5f * dt * k2[3],
                           gx, gy, gz, k3);
        
        // k4
        quaternionDerivative(qa + dt * k3[0],
                           qb + dt * k3[1],
                           qc + dt * k3[2],
                           qd + dt * k3[3],
                           gx, gy, gz, k4);
        
        // Update quaternion
        q0 += dt / 6.0f * (k1[0] + 2.0f * k2[0] + 2.0f * k3[0] + k4[0]);
        q1 += dt / 6.0f * (k1[1] + 2.0f * k2[1] + 2.0f * k3[1] + k4[1]);
        q2 += dt / 6.0f * (k1[2] + 2.0f * k2[2] + 2.0f * k3[2] + k4[2]);
        q3 += dt / 6.0f * (k1[3] + 2.0f * k2[3] + 2.0f * k3[3] + k4[3]);
        
        // Normalize quaternion
        recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
        q0 *= recipNorm;
        q1 *= recipNorm;
        q2 *= recipNorm;
        q3 *= recipNorm;
    }
}

//---------------------------------------------------------------------------------------------------
// Quaternion derivative calculation

void quaternionDerivative(float q0, float q1, float q2, float q3,
                         float wx, float wy, float wz, float* dq) {
    // Quaternion derivative from angular velocity
    dq[0] = 0.5f * (-q1 * wx - q2 * wy - q3 * wz);
    dq[1] = 0.5f * ( q0 * wx + q2 * wz - q3 * wy);
    dq[2] = 0.5f * ( q0 * wy - q1 * wz + q3 * wx);
    dq[3] = 0.5f * ( q0 * wz + q1 * wy - q2 * wx);
}

//---------------------------------------------------------------------------------------------------
// Get Euler angles from quaternion

void EnhancedMadgwickGetEulerAngles(float* roll, float* pitch, float* yaw) {
    // Convert quaternion to Euler angles (in radians)
    *roll = atan2f(2.0f * (q0 * q1 + q2 * q3), 1.0f - 2.0f * (q1 * q1 + q2 * q2));
    *pitch = asinf(2.0f * (q0 * q2 - q3 * q1));
    *yaw = atan2f(2.0f * (q0 * q3 + q1 * q2), 1.0f - 2.0f * (q2 * q2 + q3 * q3));
    
    // Convert to degrees if needed
    *roll *= 180.0f / M_PI;
    *pitch *= 180.0f / M_PI;
    *yaw *= 180.0f / M_PI;
}

//---------------------------------------------------------------------------------------------------
// Reset quaternion to initial state

void EnhancedMadgwickReset(void) {
    q0 = 1.0f;
    q1 = 0.0f;
    q2 = 0.0f;
    q3 = 0.0f;
    integralFBx = 0.0f;
    integralFBy = 0.0f;
    integralFBz = 0.0f;
}

//---------------------------------------------------------------------------------------------------
// Set algorithm parameters

void EnhancedMadgwickSetBeta(float newBeta) {
    beta = newBeta;
}

void EnhancedMadgwickSetPIGains(float newKp, float newKi) {
    Kp = newKp;
    Ki = newKi;
}

//---------------------------------------------------------------------------------------------------
// Fast inverse square-root

float invSqrt(float x) {
    float halfx = 0.5f * x;
    float y = x;
    long i = *(long*)&y;
    i = 0x5f3759df - (i>>1);
    y = *(float*)&i;
    y = y * (1.5f - (halfx * y * y));
    return y;
}

//====================================================================================================
// END OF CODE
//====================================================================================================