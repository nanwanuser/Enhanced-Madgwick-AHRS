# Enhanced Madgwick AHRS

Madgwick's AHRS algorithm with 4th-order Runge-Kutta integration for improved numerical accuracy.

## Key Improvements

- **Original Algorithm**: Gradient descent sensor fusion (unchanged)
- **Enhanced Integration**: 4th-order Runge-Kutta vs 1st-order Euler
- **Better Accuracy**: O(dt⁵) vs O(dt²) numerical error
- **Lower Sample Rate**: Works well at >50Hz (vs >400Hz recommended for original)

## Quick Start

```c
#include "EnhancedMadgwickAHRS.h"

// Update orientation (512Hz default)
EnhancedMadgwickAHRSupdate(gx, gy, gz,    // gyro (rad/s)
                           ax, ay, az,    // accel (any unit)
                           mx, my, mz);   // mag (any unit)

// Extract Euler angles
float roll  = atan2f(2.0f * (q0 * q1 + q2 * q3), 1.0f - 2.0f * (q1 * q1 + q2 * q2));
float pitch = asinf(2.0f * (q0 * q2 - q3 * q1));
float yaw   = atan2f(2.0f * (q0 * q3 + q1 * q2), 1.0f - 2.0f * (q2 * q2 + q3 * q3));
```

## Configuration

```c
// Adjust convergence rate (default: 0.1)
extern volatile float beta;
beta = 0.033f;  // Slower, more stable
beta = 0.2f;    // Faster, more responsive

// Change sample rate in source file
#define sampleFreq 100.0f  // Hz
```

## API

```c
// 9-axis fusion (gyro + accel + mag)
void EnhancedMadgwickAHRSupdate(float gx, float gy, float gz, 
                                float ax, float ay, float az, 
                                float mx, float my, float mz);

// 6-axis fusion (gyro + accel)
void EnhancedMadgwickAHRSupdateIMU(float gx, float gy, float gz, 
                                   float ax, float ay, float az);

// Global quaternion
extern volatile float q0, q1, q2, q3;
```

## Performance

| Metric | Original | Enhanced |
|--------|----------|----------|
| Integration | Euler | RK4 |
| Min Sample Rate | 400Hz | 50Hz |
| Computation | 1× | ~4× |
| Accuracy | Good | Excellent |

## License

MIT License

---

Enhanced by Nanwan (2025) | Based on S.O.H. Madgwick (2011)

---

⭐ If this enhanced implementation helps your project, please consider giving it a star!