# Global Adaptive Generative Adjustment for Generalized Linear Models

GAGA is a parameter estimation method, and it has efficient model selection ability. Its hyperparameters are easy to set. In general, alpha can be set to 1, 2 or 3, but when the collinearity of the load matrix is serious, the hyperparameters can be selected larger, such as 5. In the comparative experiment of dealing with linear models, the GAGA algorithm is competitive with adaptive lasso, SCAD, MCP. In the comparative experiment of dealing with generalized linear models, the GAGA algorithm is competitive with glmnet.

At present, GAGA can deal with linear regression, logistic and multinomial regression models, Poisson regression, Cox model.

For more detailed information, see Bin Wang, Xiaofei Wang and Jianhua Guo (2022) \<arXiv:1911.00658\>. This paper provides the theoretical properties of Gaga linear model when the load matrix is orthogonal. Further study is going on for the nonorthogonal cases and generalized linear models. These works are in part supported by the National Natural Foundation of China (No.12171076).
