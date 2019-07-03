import numpy as np
import matplotlib.pyplot as plt

style = "heart" # circle, heart
fill = False
hight = 3000

cen_lon = -45  # center longitude
cen_lat =  40  # center latitude
size = 30      # max size in deg
resolution = 300

if style == "circle":
	theta = np.linspace(0*np.pi/180., 360*np.pi/180., resolution, dtype=np.float64) 
	
	if(fill):
		data = np.full( ( resolution*resolution,3 ), -999, dtype=np.float32 )
		sizes = np.linspace(size/2, 0, resolution, dtype=np.float64) 
		iii = 0
		for ii in range(resolution):
			x = sizes[ii]*np.sin(theta) + cen_lon
			y = sizes[ii]*np.cos(theta) + cen_lat
			data[iii:iii+sizes.shape[0],0] = x
			data[iii:iii+sizes.shape[0],1] = y
			data[iii:iii+sizes.shape[0],2] = hight
			iii = iii + sizes.shape[0]

	else:
		data = np.full( ( resolution,3 ), -999, dtype=np.float32 )

		x = size/2.*np.sin(theta) + cen_lon
		y = size/2.*np.cos(theta) + cen_lat

		data[:,0] = x
		data[:,1] = y
		data[:,2] = hight

elif style == "heart":
	theta = np.linspace(0*np.pi/180., 360*np.pi/180., resolution, dtype=np.float64) 
	
	if(fill):
		data = np.full( ( resolution*resolution,3 ), -999, dtype=np.float32 )
		sizes = np.linspace(size/2, 0, resolution, dtype=np.float64) 
		iii = 0
		for ii in range(resolution):
			r = sizes[ii]
			x = r * np.sin(theta)**3 + cen_lon
			y = r * np.cos(theta)-r/2.*np.cos(2*theta)-r/6*np.cos(3*theta)-r/10.*np.cos(4*theta) + cen_lat
			data[iii:iii+sizes.shape[0],0] = x
			data[iii:iii+sizes.shape[0],1] = y
			data[iii:iii+sizes.shape[0],2] = hight
			iii = iii + sizes.shape[0]

	else:
		data = np.full( ( resolution,3 ), -999, dtype=np.float32 )

		r = size/2
		x = r * np.sin(theta)**3 + cen_lon
		r = r/2
		y = r * np.cos(theta)-r/2.*np.cos(2*theta)-r/6*np.cos(3*theta)-r/10.*np.cos(4*theta) + cen_lat

		data[:,0] = x
		data[:,1] = y
		data[:,2] = hight



print(data.shape)
np.savetxt('init_particles_generator.dat', data, fmt='%.2f %.2f %.2f',delimiter=' ',header="Generated in init_particles_generator.py")

fig = plt.figure()
plt.scatter(data[:,0], data[:,1], color='black', s=[10.])
# plt.plot(data[:,0], data[:,1])
plt.show()