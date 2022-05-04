"""
tip speed ratio = (rotaional velocity * tip distance)/wind speed

chord length = length of cross section

lift = perpendicular to the wind
drag = with the wind direction

average radius of turbine = 5-10 meters
average chord = 0.75 - 0.3 meters
omega =  72 * pi/30 radians/sec
theta = 20 to -0.1 degrees
Blades = 2 or 3

# exec '/Users/westonpace/Desktop/Julia-1.6.app/Contents/Resources/julia/bin/julia'
using Plots
using DataFrames
import Roots
using Roots
using CCBlade 


struct values
   Vinf :: Float64
   r :: Float64
   Omega :: Float64
   theta :: Float64
   Rtip :: Float64
   B :: Float64
   c :: Float64
end

constants = values(7.0 , 3.0 , 72 * pi/30 , 10 * pi/180 , 5.029 , 2 , 0.3)

function Residual(phi)
   sphi = sin(phi)
   cphi = cos(phi)
   
   alpha = constants.theta - phi
   # Re = rho * V * D / mu
   # cl= NACA12(Re,alpha)
   cl = 2*alpha
   cd = 1 + alpha^2
   cn = cl*cphi - cd*sphi
   ct = cl*sphi +cd*cphi
   # Tip Correction
   smf = constants.B/2*(constants.Rtip-constants.r)/(constants.r*abs(sphi))
   F = 2/pi * acos(exp(-smf))
   lam = constants.Omega * constants.r /constants.Vinf
   sigp = constants.B*constants.c/(2*pi*constants.r)
   a = 1/((4*F*sphi^2/(sigp*cn))-1)
   ap = 1/((4*F*sphi*cphi/(sigp*ct))+1)
   R = sphi/(1+a) - cphi/(lam*(1-ap))
   return R
end


f(phi)= Residual(phi)
a=Roots.find_zero(f,(0+1e-7,pi/2))
print(a * 180/pi)

"""
#Practice BEM
"""
function Residual(phi, Vinf, B, c, r, Rtip, Omega, Theta)
   sphi = sin(phi)
   cphi = cos(phi)
   alpha = theta - phi
   cl = 2*alpha
   cd = 1 + alpha^2
   cn = cl*cphi - cd*sphi
   ct = sl*sphi + cdcphi
   smf = B/2*(Rtip-r)/(r*abs(sphi))
   F = 2/pi * acos(exp(-smf))
   lam = Omega * r / Vinf
   sigp = B*c/(2*pi*r)
   a = 1/((4*F*sphi^2)/(sigp*cn))-1)
   ap = 1/((4*F*sphi*cphi)/(sigp*ct))+1)
   R = (sphi/(1+a))-(cphi/(lam(1-ap)))
   return R
end
"""
"""

#for power curve with 5 section
#seet up residual function and try to change inputs during each loop
using Plots
using DataFrames
using Roots
radius = [2.0,4.0,6.0,8.0,10.0]
chord = [0.75, 0.637 , 0.525 , 0.412 , 0.3]
theta1 = [20.0 , 14.975 , 9.95 , 4.925 , -0.1]
container = []
Final = []
X = []
function Residual(phi , theta2 , r , c , b , w , Rtip , Vinf)
   sphi = sin(phi)
   cphi = cos(phi)
   
   alpha = theta2*(pi/180) - phi
   # Re = rho * V * D / mu
   # cl= NACA12(Re,alpha)
   cl = 2*alpha
   cd = 1 + alpha^2
   cn = cl*cphi - cd*sphi
   ct = cl*sphi +cd*cphi
   # Tip Correction
   smf = b/2*(Rtip-r)/(r*abs(sphi))
   F = 2/pi * acos(exp(-smf))
   lam = w * r /Vinf
   sigp = b*c/(2*pi*r)
   a = 1/((4*F*sphi^2/(sigp*cn))-1)
   ap = 1/((4*F*sphi*cphi/(sigp*ct))+1)
   R = sphi/(1+a) - cphi/(lam*(1-ap))
   return R
end

function Solved_Residual(a, theta2 , r , c , b , w , Rtip , Vinf)
   sphi = sin(a)
   cphi = cos(a)
   
   alpha = theta2*(pi/180) - a
   # Re = rho * V * D / mu
   # cl= NACA12(Re,alpha)
   cl = 2*alpha
   cd = 1 + alpha^2
   cn = cl*cphi - cd*sphi
   ct = cl*sphi +cd*cphi
   # Tip Correction
   smf = b/2*(Rtip-r)/(r*abs(sphi))
   F = 2/pi * acos(exp(-smf))
   lam = w * r /Vinf
   sigp = b*c/(2*pi*r)
   a = 1/((4*F*sphi^2/(sigp*cn))-1)
   ap = 1/((4*F*sphi*cphi/(sigp*ct))+1)
   R = sphi/(1+a) - cphi/(lam*(1-ap))
   return ap
end

for Vinf in 7.0:20.0
   for w in 1:5
      r = radius[w]
      c = chord[w]
      theta2 = theta1[w]
      b = 2.0
      w = 50 * pi/30
      Rtip =  10.1
      f(phi) = Residual(phi, theta2 , r , c , b , w , Rtip , Vinf)
      a = Roots.find_zero(f , (0+1e-7,pi/2))
      ct = 2*(theta2*(pi/180) - a)*sin(a) + (1+(theta2*(pi/180) - a)^2)*cos(a)
      Q = b*ct*(1/2)*1.225*((Vinf*(1+a))^2+(w*r*(1-Solved_Residual(a, theta2 , r , c , b , w , Rtip , Vinf)))^2)*c*r
      Q2 = Q*2
      P = Q2 * w
      push!(container , P)
   end
   P_Total = sum(container)
   push!(Final , P_Total)
   empty!(container)
   push!(X , Vinf)
end
"""



using Plots
using DataFrames
using Roots

alpha = [-20.1000 ,-18.1000 ,-16.1000 ,-14.2000 ,-12.2000 ,-10.1000 , -8.2000 , -6.1000 
, -4.1000 ,-2.1000 , 0.1000 ,  2.0000 , 4.1000 ,6.2000  ,8.1000 , 10.2000 , 
11.3000  , 12.1000 , 13.2000 , 14.2000  ,15.3000     , 16.3000    , 17.1000     , 18.1000     
, 19.1000     , 20.1000    , 22.0000     ,24.1000      ,26.2000]

cl =[-0.56,-0.67,-0.79,-0.84,-0.7,-0.63,-0.56,-0.64,-0.42,-0.21,0.05,0.3,0.54,0.79,0.9,0.93,0.92,
0.95,0.99,1.01,1.02,1,0.94,0.85,0.7,0.66,0.7,0.79,0.88]

cd = [0.3027,0.3069,0.1928,0.0898,0.0553,0.039,0.0872,0.0231,0.0117,0.0107,0.0142,0.0163,0.0188,0.0216,
0.0266,0.071,0.0303,0.0369,0.0509,0.0648,0.0776,0.0917,0.0994,0.2306,0.3142,0.3186,0.3694,0.4457,
0.526]

chord = []
for i in 0:28
   x = 0.75+i*(0.3-0.75)/28
   push!(chord , x)
end

radius= []
for i in 0:28
   x = 0.1+i*(5.029-0.1)/28
   push!(radius , x)
end
container = []
Final = []
X = []
function Residual(phi , theta2 , r , c , b , Rtip , Vinf , lift, drag)
   sphi = sin(phi)
   cphi = cos(phi)
   
   alpha = theta2*(pi/180) - phi
   # Re = rho * V * D / mu
   # cl= NACA12(Re,alpha)
   cl = lift
   cd = drag
   cn = cl*cphi - cd*sphi
   ct = cl*sphi +cd*cphi
   # Tip Correction
   smf = b/2*(Rtip-r)/(r*abs(sphi))
   F = 2/pi * acos(exp(-smf))
   lam = 7
   sigp = b*c/(2*pi*r)
   a = 1/((4*F*sphi^2/(sigp*cn))-1)
   ap = 1/((4*F*sphi*cphi/(sigp*ct))+1)
   R = sphi/(1+a) - cphi/(lam*(1-ap))
   return R
end

function Solved_Residual_1(a, theta2 , r , c , b , Rtip , Vinf, lift, drag)
   sphi = sin(a)
   cphi = cos(a)
   
   alpha = theta2*(pi/180) - a
   # Re = rho * V * D / mu
   # cl= NACA12(Re,alpha)
   cl = lift
   cd = drag
   cn = cl*cphi - cd*sphi
   ct = cl*sphi +cd*cphi
   # Tip Correction
   smf = b/2*(Rtip-r)/(r*abs(sphi))
   F = 2/pi * acos(exp(-smf))
   lam = 7
   w = 7*Vinf/r
   sigp = b*c/(2*pi*r)
   a = 1/((4*F*sphi^2/(sigp*cn))-1)
   ap = 1/((4*F*sphi*cphi/(sigp*ct))+1)
   R = sphi/(1+a) - cphi/(lam*(1-ap))
   return ap
end

function Solved_Residual_2(a, theta2 , r , c , b , Rtip , Vinf, lift, drag)
   sphi = sin(a)
   cphi = cos(a)
   
   alpha = theta2*(pi/180) - a
   # Re = rho * V * D / mu
   # cl= NACA12(Re,alpha)
   cl = lift
   cd = drag
   cn = cl*cphi - cd*sphi
   ct = cl*sphi +cd*cphi
   # Tip Correction
   smf = b/2*(Rtip-r)/(r*abs(sphi))
   F = 2/pi * acos(exp(-smf))
   lam = 7
   w = 7*Vinf/r
   sigp = b*c/(2*pi*r)
   a = 1/((4*F*sphi^2/(sigp*cn))-1)
   ap = 1/((4*F*sphi*cphi/(sigp*ct))+1)
   R = sphi/(1+a) - cphi/(lam*(1-ap))
   return w
end

for Vinf in 1.0:0.01:14.0
   for i in 1:24
      lift = cl[i]
      drag = cd[i]
      r = radius[i]
      c = chord[i]
      theta2 = alpha[i]
      b = 2.0
      Rtip =  5.03
      f(phi) = Residual(phi, theta2 , r , c , b , Rtip , Vinf, lift, drag)
      a = Roots.find_zero(f , (-pi/2,pi/2))
      ct = 2*(theta2*(pi/180) - a)*sin(a) + (1+(theta2*(pi/180) - a)^2)*cos(a)
      Q = b*ct*(1/2)*1.225*((Vinf*(1+a))^2+(Solved_Residual_2(a, theta2 , r , c , b , Rtip , Vinf, lift, drag)*r*(1-Solved_Residual_1(a, theta2 , r , c , b , Rtip , Vinf, lift, drag)))^2)*c*r
      Q2 = Q*0.176
      P = Q2 * Solved_Residual_2(a, theta2 , r , c , b , Rtip , Vinf, lift, drag)
      push!(container , P)
   end
   P_Total = sum(container)
   push!(Final , P_Total/1000)
   empty!(container)
   push!(X , Vinf)
end


Power =[]
axis =[]
for i in 1:0.01:14
   P = 0.4*(1/2)*(pi*5.03^2)*1.225*i^3
   push!(Power, P/1000)
   push!(axis,i)
end

Percent = []
Place_Holder = []
for i in 1:length(Final)
   push!(Percent, (Final[i]-Power[i]))
   push!(Place_Holder, (Final[i]-Power[i])*0.01)
end
Average = sum(Place_Holder)*(1/(20-7))

for i in 1:1301
   push!(thing , Average)
end