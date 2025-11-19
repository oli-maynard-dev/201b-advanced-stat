import math
import numpy as np

# Data
n = 5146  
x = 1717  
m = 6099  
y = 2030  

# Log f(x,y|H_1)
log_f_H1 = -np.log(n + 1) - np.log(m + 1)

# Log f(x,y|H_0)
log_binom_n_x = math.lgamma(n + 1) - math.lgamma(x + 1) - math.lgamma(n - x + 1)
log_binom_m_y = math.lgamma(m + 1) - math.lgamma(y + 1) - math.lgamma(m - y + 1)
log_beta = (math.lgamma(x + y + 1) + 
            math.lgamma(n + m - x - y + 1) - 
            math.lgamma(n + m + 2))

log_f_H0 = log_binom_n_x + log_binom_m_y + log_beta

# Bayes Factor
log_BF10 = log_f_H1 - log_f_H0
BF10 = np.exp(log_BF10)

# Results
print(f"log(BF_10) = {log_BF10:.6f}")
print(f"BF_10 = {BF10:.6f}")