//=====================================================================================================
// EnhancedMadgwickAHRS.c
//=====================================================================================================
//
// Enhanced implementation of Madgwick's IMU and AHRS algorithms.
// Based on the original Madgwick algorithm with 4th-order Runge-Kutta integration.
// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
//
// Date			Author          Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
// 19/02/2012	SOH Madgwick	Magnetometer measurement is normalised
// 04/07/2025	Nanwan          Added 4th-order Runge-Kutta integration
//
//=====================================================================================================

//---------------------------------------------------------------------------------------------------
// Header files

#include "EnhancedMadgwickAHRS.h"
#include <math.h>

//---------------------------------------------------------------------------------------------------
// Definitions

#define sampleFreq	512.0f		// sample frequency in Hz
#define betaDef		0.1f		// 2 * proportional gain

//---------------------------------------------------------------------------------------------------
// Variable definitions

volatile float beta = betaDef;								// 2 * proportional gain (Kp)
volatile float q0 = 1.0f, q1 = 0.0f, q2 = 0.0f, q3 = 0.0f;	// quaternion of sensor frame relative to auxiliary frame

//---------------------------------------------------------------------------------------------------
// Function declarations

float invSqrt(float x);

//====================================================================================================
// Functions

//---------------------------------------------------------------------------------------------------
// AHRS algorithm update

void EnhancedMadgwickAHRSupdate(float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz) {
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float hx, hy;
	float _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;
	
	// Runge-Kutta variables
	float k1[4], k2[4], k3[4], k4[4];
	float q_temp[4];
	float dt = 1.0f / sampleFreq;
	float gx_corrected, gy_corrected, gz_corrected;

	// Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
	if((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f)) {
		EnhancedMadgwickAHRSupdateIMU(gx, gy, gz, ax, ay, az);
		return;
	}

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
	qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
	qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
	qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;   

		// Normalise magnetometer measurement
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
		_2bx = sqrt(hx * hx + hy * hy);
		_2bz = -_2q0mx * q2 + _2q0my * q1 + mz * q0q0 + _2q1mx * q3 - mz * q1q1 + _2q2 * my * q3 - mz * q2q2 + mz * q3q3;
		_4bx = 2.0f * _2bx;
		_4bz = 2.0f * _2bz;

		// Gradient decent algorithm corrective step
		s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q3 + _2bz * q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q2 + _2bz * q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q3 - _4bz * q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * q2 - _2bz * q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * q1 + _2bz * q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * q0 - _4bz * q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * q3 + _2bz * q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * q0 + _2bz * q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= beta * s0;
		qDot2 -= beta * s1;
		qDot3 -= beta * s2;
		qDot4 -= beta * s3;
		
		// Store corrected gyroscope values for RK4
		gx_corrected = gx - 2.0f * beta * s1;
		gy_corrected = gy - 2.0f * beta * s2;
		gz_corrected = gz - 2.0f * beta * s3;
	} else {
		gx_corrected = gx;
		gy_corrected = gy;
		gz_corrected = gz;
	}

	// 4th-order Runge-Kutta integration
	// k1 = f(t, y)
	k1[0] = qDot1;
	k1[1] = qDot2;
	k1[2] = qDot3;
	k1[3] = qDot4;
	
	// k2 = f(t + dt/2, y + dt/2 * k1)
	q_temp[0] = q0 + 0.5f * dt * k1[0];
	q_temp[1] = q1 + 0.5f * dt * k1[1];
	q_temp[2] = q2 + 0.5f * dt * k1[2];
	q_temp[3] = q3 + 0.5f * dt * k1[3];
	
	k2[0] = 0.5f * (-q_temp[1] * gx_corrected - q_temp[2] * gy_corrected - q_temp[3] * gz_corrected);
	k2[1] = 0.5f * ( q_temp[0] * gx_corrected + q_temp[2] * gz_corrected - q_temp[3] * gy_corrected);
	k2[2] = 0.5f * ( q_temp[0] * gy_corrected - q_temp[1] * gz_corrected + q_temp[3] * gx_corrected);
	k2[3] = 0.5f * ( q_temp[0] * gz_corrected + q_temp[1] * gy_corrected - q_temp[2] * gx_corrected);
	
	// k3 = f(t + dt/2, y + dt/2 * k2)
	q_temp[0] = q0 + 0.5f * dt * k2[0];
	q_temp[1] = q1 + 0.5f * dt * k2[1];
	q_temp[2] = q2 + 0.5f * dt * k2[2];
	q_temp[3] = q3 + 0.5f * dt * k2[3];
	
	k3[0] = 0.5f * (-q_temp[1] * gx_corrected - q_temp[2] * gy_corrected - q_temp[3] * gz_corrected);
	k3[1] = 0.5f * ( q_temp[0] * gx_corrected + q_temp[2] * gz_corrected - q_temp[3] * gy_corrected);
	k3[2] = 0.5f * ( q_temp[0] * gy_corrected - q_temp[1] * gz_corrected + q_temp[3] * gx_corrected);
	k3[3] = 0.5f * ( q_temp[0] * gz_corrected + q_temp[1] * gy_corrected - q_temp[2] * gx_corrected);
	
	// k4 = f(t + dt, y + dt * k3)
	q_temp[0] = q0 + dt * k3[0];
	q_temp[1] = q1 + dt * k3[1];
	q_temp[2] = q2 + dt * k3[2];
	q_temp[3] = q3 + dt * k3[3];
	
	k4[0] = 0.5f * (-q_temp[1] * gx_corrected - q_temp[2] * gy_corrected - q_temp[3] * gz_corrected);
	k4[1] = 0.5f * ( q_temp[0] * gx_corrected + q_temp[2] * gz_corrected - q_temp[3] * gy_corrected);
	k4[2] = 0.5f * ( q_temp[0] * gy_corrected - q_temp[1] * gz_corrected + q_temp[3] * gx_corrected);
	k4[3] = 0.5f * ( q_temp[0] * gz_corrected + q_temp[1] * gy_corrected - q_temp[2] * gx_corrected);
	
	// Update quaternion: y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
	q0 += dt / 6.0f * (k1[0] + 2.0f * k2[0] + 2.0f * k3[0] + k4[0]);
	q1 += dt / 6.0f * (k1[1] + 2.0f * k2[1] + 2.0f * k3[1] + k4[1]);
	q2 += dt / 6.0f * (k1[2] + 2.0f * k2[2] + 2.0f * k3[2] + k4[2]);
	q3 += dt / 6.0f * (k1[3] + 2.0f * k2[3] + 2.0f * k3[3] + k4[3]);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
}

//---------------------------------------------------------------------------------------------------
// IMU algorithm update

void EnhancedMadgwickAHRSupdateIMU(float gx, float gy, float gz, float ax, float ay, float az) {
	float recipNorm;
	float s0, s1, s2, s3;
	float qDot1, qDot2, qDot3, qDot4;
	float _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;
	
	// Runge-Kutta variables
	float k1[4], k2[4], k3[4], k4[4];
	float q_temp[4];
	float dt = 1.0f / sampleFreq;
	float gx_corrected, gy_corrected, gz_corrected;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-q1 * gx - q2 * gy - q3 * gz);
	qDot2 = 0.5f * (q0 * gx + q2 * gz - q3 * gy);
	qDot3 = 0.5f * (q0 * gy - q1 * gz + q3 * gx);
	qDot4 = 0.5f * (q0 * gz + q1 * gy - q2 * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = invSqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;   

		// Auxiliary variables to avoid repeated arithmetic
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

		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		s2 = 4.0f * q0q0 * q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		s3 = 4.0f * q1q1 * q3 - _2q1 * ax + 4.0f * q2q2 * q3 - _2q2 * ay;
		recipNorm = invSqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= beta * s0;
		qDot2 -= beta * s1;
		qDot3 -= beta * s2;
		qDot4 -= beta * s3;
		
		// Store corrected gyroscope values for RK4
		gx_corrected = gx - 2.0f * beta * s1;
		gy_corrected = gy - 2.0f * beta * s2;
		gz_corrected = gz - 2.0f * beta * s3;
	} else {
		gx_corrected = gx;
		gy_corrected = gy;
		gz_corrected = gz;
	}

	// 4th-order Runge-Kutta integration
	// k1 = f(t, y)
	k1[0] = qDot1;
	k1[1] = qDot2;
	k1[2] = qDot3;
	k1[3] = qDot4;
	
	// k2 = f(t + dt/2, y + dt/2 * k1)
	q_temp[0] = q0 + 0.5f * dt * k1[0];
	q_temp[1] = q1 + 0.5f * dt * k1[1];
	q_temp[2] = q2 + 0.5f * dt * k1[2];
	q_temp[3] = q3 + 0.5f * dt * k1[3];
	
	k2[0] = 0.5f * (-q_temp[1] * gx_corrected - q_temp[2] * gy_corrected - q_temp[3] * gz_corrected);
	k2[1] = 0.5f * ( q_temp[0] * gx_corrected + q_temp[2] * gz_corrected - q_temp[3] * gy_corrected);
	k2[2] = 0.5f * ( q_temp[0] * gy_corrected - q_temp[1] * gz_corrected + q_temp[3] * gx_corrected);
	k2[3] = 0.5f * ( q_temp[0] * gz_corrected + q_temp[1] * gy_corrected - q_temp[2] * gx_corrected);
	
	// k3 = f(t + dt/2, y + dt/2 * k2)
	q_temp[0] = q0 + 0.5f * dt * k2[0];
	q_temp[1] = q1 + 0.5f * dt * k2[1];
	q_temp[2] = q2 + 0.5f * dt * k2[2];
	q_temp[3] = q3 + 0.5f * dt * k2[3];
	
	k3[0] = 0.5f * (-q_temp[1] * gx_corrected - q_temp[2] * gy_corrected - q_temp[3] * gz_corrected);
	k3[1] = 0.5f * ( q_temp[0] * gx_corrected + q_temp[2] * gz_corrected - q_temp[3] * gy_corrected);
	k3[2] = 0.5f * ( q_temp[0] * gy_corrected - q_temp[1] * gz_corrected + q_temp[3] * gx_corrected);
	k3[3] = 0.5f * ( q_temp[0] * gz_corrected + q_temp[1] * gy_corrected - q_temp[2] * gx_corrected);
	
	// k4 = f(t + dt, y + dt * k3)
	q_temp[0] = q0 + dt * k3[0];
	q_temp[1] = q1 + dt * k3[1];
	q_temp[2] = q2 + dt * k3[2];
	q_temp[3] = q3 + dt * k3[3];
	
	k4[0] = 0.5f * (-q_temp[1] * gx_corrected - q_temp[2] * gy_corrected - q_temp[3] * gz_corrected);
	k4[1] = 0.5f * ( q_temp[0] * gx_corrected + q_temp[2] * gz_corrected - q_temp[3] * gy_corrected);
	k4[2] = 0.5f * ( q_temp[0] * gy_corrected - q_temp[1] * gz_corrected + q_temp[3] * gx_corrected);
	k4[3] = 0.5f * ( q_temp[0] * gz_corrected + q_temp[1] * gy_corrected - q_temp[2] * gx_corrected);
	
	// Update quaternion: y_new = y + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
	q0 += dt / 6.0f * (k1[0] + 2.0f * k2[0] + 2.0f * k3[0] + k4[0]);
	q1 += dt / 6.0f * (k1[1] + 2.0f * k2[1] + 2.0f * k3[1] + k4[1]);
	q2 += dt / 6.0f * (k1[2] + 2.0f * k2[2] + 2.0f * k3[2] + k4[2]);
	q3 += dt / 6.0f * (k1[3] + 2.0f * k2[3] + 2.0f * k3[3] + k4[3]);

	// Normalise quaternion
	recipNorm = invSqrt(q0 * q0 + q1 * q1 + q2 * q2 + q3 * q3);
	q0 *= recipNorm;
	q1 *= recipNorm;
	q2 *= recipNorm;
	q3 *= recipNorm;
}

//---------------------------------------------------------------------------------------------------
// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root

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