# Used to test ideas
import numpy as np
def near(n):
    return np.ceil(n - 0.5) + 0.5

print(near(2.50001))