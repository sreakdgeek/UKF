# Unscented Kalman Filters

Extended Kalman filters were quite useful to track an object that typically has linear motion. In EKF, object's state distribution is approximated by a Gaussian
Random Variable. When such a Gaussian random variable is mapped to a non-linear function, then the result is no longer a Gaussian distribution. Hence Extended
Kalman Filters use Taylor series to linearize the non-linear function. It has been shown that this can introduce large errors in the posterior mean and covariance
of the transformed Gaussian random variable. Unscented Kalman Filters overcome this limitation by choosing sample points that completely capture the true mean and
covariance of the Gaussian random variable which are accurate up to 3rd order of Taylor series exapansion. More details of the idea and implementation can be found
in the original paper.

[//]: # (Image References) 
[image1]: ./images/UKF_Roadmap.png
[image2]: ./images/UKF_Viz_DS1.png
[image3]: ./images/UKF_Viz_DS2.png
[image4]: ./images/Results.png

Original paper is described here: https://www.seas.harvard.edu/courses/cs281/papers/unscented.pdf

### UKF Road Map

[image1]

---

### Results

#### Graph for Dataset 1

[image2]


#### Graph for Dataset 2

[image3]

---

#### RMSE

[image4]

