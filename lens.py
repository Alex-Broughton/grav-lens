import matplotlib.pyplot as plt
import numpy as np
import scipy.integrate as integrate
from scipy.optimize import fsolve

# Constants
Om = 0.3
Ol = 0.7
H0 = 70.0 # km/s/Mpc
c = 3.0e8 # m/s
G = 6.67e-11 # kg m^3 / s^2

z_source = 1.66 # Input parameters
z_image = 0.41

image1 = [0.237,0.798]
image2 = [0.005,-0.610]
source = [0.0,0.0]

plt.figure(1) # Find Einstein radius
plt.figure(figsize=(4,4))
plt.scatter([image1[0]], [image1[1]], c="#000000", s=20)
plt.scatter([image2[0]], [image2[1]], c="#000000", s=20)
einsteinCircle = plt.Circle((0.0, 0.0), 0.720, color='red', alpha = 0.2)
ax=plt.gca()
ax.add_patch(einsteinCircle)
plt.scatter([source[0]], [source[1]], c="#000000", s=60)
plt.xlim(-1.0,1.0)
plt.ylim(-1.0,1.0)
plt.xlabel("dec [arcseconds]")
plt.ylabel("az [arcseconds]")
print("Given distribution:")
plt.show()

MpcToM = 3.08567758e22
arcsecToRad = 4.8481e-6
kgToSM = 1.0 / 1.9884e30
einsteinRadius_arc = 0.720 # arcseconds
einsteinRadius_rad = einsteinRadius_arc * arcsecToRad

# From here, we can write the gravitational lens equation:
# theta = beta + alpha

def E(z): # Calculate E(z)
    E = np.sqrt((Om * (1.0 + z)**3) + Ol)
    return E

def integrand(z): # Comoving distance
    return 3000. / ( (H0/100.) * E(z))

def Dc(z): # Comoving distance
    sol = integrate.quad(integrand, 0.0, z)
    return sol[0]

def getDds(zs, zd): #Dds
    Ds = Dc(zs)
    Dd = Dc(zd)
    Dds = (Ds - Dd) / (1 + zs)
    return Dds

def getD(z): # Distance in expanding universe
    return Dc(z) / (1 + z)

Dds = getDds(z_source, z_image) * MpcToM # Convert to meters
Dd = getD(z_image) * MpcToM
Ds = getD(z_source) * MpcToM

def M(Dd, Ds, Dds, theta_E): # Calculate mass
    return ( (theta_E**2.0) * (c*c*Dd*Ds) ) / (4.0 * G * Dds)

lensMass = M(Dd, Ds, Dds, einsteinRadius_rad)
lensMass_SM = lensMass * kgToSM
print("Object mass (kg): ")
print(lensMass_SM)

# Calculate image positions

def E_radius_2(mass): # Get Einstein radius
    return (4*G*mass*Dds) / (c*c*Dd*Ds)

def func(theta,*params): # Minimization function
    beta_x, beta_y, M = params[0]
    f = [0,0]
    theta_E_2 = E_radius_2(M)
    f[0] = (theta[0] - theta_E_2*(theta[0]/(theta[0]**2 + theta[1]**2)) - beta_x)
    f[1] = (theta[1] - theta_E_2*(theta[1]/(theta[0]**2 + theta[1]**2)) - beta_y)
    
    return f



beta_x_init = 0.04 * arcsecToRad # Set initial conditions for source position .04
beta_y_init = 0.2 * arcsecToRad # .2
params = [beta_x_init, beta_y_init, lensMass] # Use lens mass
image1_rad = [0,0]
image2_rad = [0,0]
image1_rad[0] = image1[0] * arcsecToRad # Get images in radians
image1_rad[1] = image1[1] * arcsecToRad
image2_rad[0] = image2[0] * arcsecToRad
image2_rad[1] = image2[1] * arcsecToRad

image1_calculated = fsolve(func, image1_rad, params) / arcsecToRad # Solve for image positions
image2_calculated = fsolve(func, image2_rad, params) / arcsecToRad
    
plt.figure(2)                  # Plot solution
plt.figure(figsize=(4,4))
plt.xlim(-1.0,1.0)
plt.ylim(-1.0,1.0)
plt.scatter([image1[0]], [image1[1]], c="#000000", s=20)
plt.scatter([image2[0]], [image2[1]], c="#000000", s=20)
plt.scatter([source[0]], [source[1]], c="#000000", s=60)
plt.scatter( [ image1_calculated[0] ] , [ image1_calculated[1] ] , c="#AA00FF", marker="*", s=50)
plt.scatter( [ image2_calculated[0] ] , [ image2_calculated[1] ] , c="#AA00FF", marker="*", s=50)
plt.scatter([beta_x_init / arcsecToRad], [beta_y_init / arcsecToRad], c="#AA00FF", marker="*", s=150)
plt.xlabel("dec [arcsec]")
plt.ylabel("az [arcsec]")
print("Calculated position: ")
plt.show()
