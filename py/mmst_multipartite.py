# mmst_multipartite.py
# Mathew Titus, Sunstrand Technical Consulting
# May, 2024
# 
# Computing combinatorial quantities related to the mean size
# of minimal spanning trees among regular k-partite graphs with
# unequal partition weights.
# 
########################################################################

import numpy as np
import matplotlib.pyplot as plt
from itertools import combinations_with_replacement as cwr
from functools import reduce
import sys
import json
# from argparse import ArgumentParser

# parser = ArgumentParser()
# parser.add_argument("p", type=float)

# define k - the number of particle types
k = 4
assert k == 4, f"Interaction strength matrix is not defined for k = {k}"

def factorial(n: int):
  assert isinstance(n, int), f"Factorial function only accepts integer inputs.\n(Received {type(n)})";
  if n < 0: return np.nan;
  return 1 if ((n==1)|(n==0)) else n*factorial(n-1);

  
# define function for each summand {x: <x | 1> = n}
def summand(x: np.ndarray, V: np.ndarray):
  xm = [float(num) for num in x - 1];
  y = np.matmul(V, x)**np.array(xm)

  # compute leading constant
  c0 = factorial(int(x.sum()))

  # compute product of components of y
  c1 = reduce(
      lambda a,b: a*b,
      y,
      1
    )
  
  # compute denominator, x! = x_1! * x_2! * ... * x_k!
  c2 = reduce(
      lambda a,b: a*b,
      map(
        lambda en: factorial(int(en)), 
        x
      ),
      1
    )

  return c0 * c1 / c2;


# define function to compute partitions
def get_partitions(n, k):
  '''
  number of ways to designate your n particles into k urns:
  place k-1 separators among/around the ordered n particles (n+1 choices)
  the particles above the highest separator are urn 1 members,
  the particles between the first & second sep are in urn 2, ...
  the placement of separators is done by choosing one of the n+1 spaces k-1 times
  with replacement and without regard to order; use the itertools function for this.
  '''
  raw_parts = list(
      cwr(np.arange(0, n+1), k-1)
    ) # weak

    # first bin:  rp[0]
    # middle bins:  np.diff(rp)
    # final bin:  n - rp[-1]

  # would be nice to perform this process using map, not listing
  out = np.empty((len(raw_parts), k));
  for ent, y in enumerate(raw_parts):
    out[ent, :] = [y[0]] + np.diff(y).tolist() + [n-y[-1]]

  return out;


# calculate the full LHS of (23)
def compute_lhs(n, k, V):
  lhs = 0;
  for part in get_partitions(n,k):
    lhs += summand(part, V)
  return lhs;


# expected RHS
def compute_rhs(n,k):
  rhs = k * (k-1)**(n-1) * n**(n-4)
  return rhs;


# define n - the number of total particles
n = 2

# define interaction strength matrix
p, q, r = [1, 1, 5]

resolution = 100

if __name__ == "__main__":
  args = sys.argv
  outfilename = f"{args[1]}_{args[2]}_{args[3]}_{args[4]}.csv"
  print(args)
  assert len(args) == 5
  p = float(args[1])
  q = float(args[2])
  r = float(args[3])
  n = int(args[4])

  p_domain = np.linspace(0.1, 50, resolution)
  q_domain = np.linspace(0.1, 50, resolution)
  [X,Y] = np.meshgrid(p_domain, q_domain)
  Z = np.zeros(X.shape)

  for eye, p in enumerate(p_domain):
    for jay, q in enumerate(q_domain):
      V = np.array([
        [0, p, q, r],
        [p, 0, r, q],
        [q, r, 0, p],
        [r, q, p, 0]
      ])

      lhs = compute_lhs(n,k,V)
      rhs = compute_rhs(n,k)
      Z[eye, jay] = lhs/rhs;

  with open(f"../output/{outfilename}", 'w') as f:
    json.dump(Z.tolist(), f)

# plt.plot(r_domain, F)


# ## Plot result
# from matplotlib import cm

# fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
# surf = ax.plot_surface(X, Y, (Z), cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)

# plt.show()



## Some output

# p,q,r  n :1     2     3       4         5     6

# 1,1,1  : 1.0 constant

# 1,1,2  : 

# 1,1,3  : 1/3   1/2   83/105  4141/3240 





# 105 = 3 * 5 * 7
# 4141 = 41 * 101
# 3240 = 2**3 * 3**4 * 5

