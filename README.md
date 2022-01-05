# Decorrelation-Score Function Tutorial
The "Decorrelation-Score Function Tutorial" implements a simulation function that performs some applications for Decorrelation-Score Function introduced by Ning, Liu, et al. (2017a) and Ning, Liu, et al. (2017b). This implementation can replicate some simple regression methods with decorrelation-Score Function.

Paper link:https://arxiv.org/abs/1412.8765. 
## Introdcution
In this study, I focus on sparse high dimensional modeling with the general decorrelated score methods proposed in Ning, Liu, et al. (2017a) and explore to which extent does this general theory applied in regression models. A research question arises about the applicability and generalization of decorrelated score methods for interest models and other available regression models. Linear regression, logistic regression, and Poisson regression are used to explore this question.
## Simulation Process
Throughout the simulation study, I first set the data generator process (DGP) of the covariates X: n = 100 independent and identical distribution samples with a multivariate Gaussian distribution Nd(0,Σ), where d = 100,200,500 and Σ is a diagonal-constant matrix with Σij = ρ^|i−j|. ρ has four potential values, namely, 0.25, 0.4, 0.6, and 0.75. The magnitude of ρ determines the strength of the collinearity of the data from DGP. For the true value of β∗, it satisfies ‖β∗‖0 = s, where βS = (1,...,1) is Dirac measure with s = 3. For the simulation process of linear regression, there is a standard Gaussian noise assumed in DGP, that is, Y = Xβ∗ + ε, where ε is a n ×1 vector following standard normal distribution. Regarding to generator of Y in logistic regression, Y is assumed to follow the binomial distribution with the probability of success on each trial.
## Results Display 
Here are some replication simulation results of the power with linear regression, logistic regression, and Poisson regression:
![result2](https://user-images.githubusercontent.com/59536847/147810193-edd96b10-2c96-4750-b17b-cc2012ffa2f7.PNG)
![result3](https://user-images.githubusercontent.com/59536847/147810217-e7dfe650-1617-433c-9f5a-51409d243eca.PNG)
![result_1](https://user-images.githubusercontent.com/59536847/147810221-bdc34582-e39b-46f0-a27c-fe18cc0e444a.PNG)
![result4](https://user-images.githubusercontent.com/59536847/147810222-1fc2e3cb-8cef-4c59-9f83-75bd12e86477.PNG)
## Reference 
Ning, Y., & Liu, H. (2017). A general theory of hypothesis tests and confidence regions for sparse high dimensional models. The Annals of Statistics, 45(1), 158-195.
