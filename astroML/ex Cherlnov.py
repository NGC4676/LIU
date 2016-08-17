from matplotlib.pylab import *

data = fetch_sdss_specgals()[:]
print data['mjd'][28],data['plate'][28],data['fiberID'][28]
def cface(ax, x6,x8,x12):
    # x1 = height  of upper face
    # x2 = overlap of lower face
    # x3 = half of vertical size of face
    # x4 = width of upper face
    # x5 = width of lower face
    # x6 = length of nose
    # x7 = vertical position of mouth
    # x8 = curvature of mouth
    # x9 = width of mouth
    # x10 = vertical position of eyes
    # x11 = separation of eyes
    # x12 = slant of eyes
    # x13 = eccentricity of eyes
    # x14 = size of eyes
    # x15 = position of pupils
    # x16 = vertical position of eyebrows
    # x17 = slant of eyebrows
    # x18 = size of eyebrows

    # transform some values so that input between 0,1 yields variety of output
    x1=0.9
    x2=0.6
    x3 = 0.5
    x4 = 1
    x5 = 1
    x7=0.7
    x8 = 5*(x8+.001)
    x9= 1.1*x8
    x10= x6/1.2
    x11 = (1.1-x6)/0.9
    x12 = 2*(x12-.5)
    x13 = x12
    x14 = 0.2*x8
    x15 = 0.1

    # top of face, in box with l=-x4, r=x4, t=x1, b=x3
    e = mpl.patches.Ellipse( (0,(x1+x3)/2), 2*x4, (x1-x3), fc='white', linewidth=2)
    ax.add_artist(e)

    # bottom of face, in box with l=-x5, r=x5, b=-x1, t=x2+x3
    e = mpl.patches.Ellipse( (0,(-x1+x2+x3)/2), 2*x5, (x1+x2+x3), fc='white', linewidth=2)
    ax.add_artist(e)

    # cover overlaps
    e = mpl.patches.Ellipse( (0,(x1+x3)/2), 2*x4, (x1-x3), fc='white', ec='none')
    ax.add_artist(e)
    e = mpl.patches.Ellipse( (0,(-x1+x2+x3)/2), 2*x5, (x1+x2+x3), fc='white', ec='none')
    ax.add_artist(e)
    
    # draw nose
    plot([0,0], [-x6/2, x6/2], 'k')
    
    # draw mouth
    p = mpl.patches.Arc( (0,-x7+.5/x8), 1/x8, 1/x8, theta1=270-180/pi*arctan(x8*x9), theta2=270+180/pi*arctan(x8*x9))
    ax.add_artist(p)
    
    # draw eyes
    p = mpl.patches.Ellipse( (-x11-x14/2,x10), x14, x13*x14, angle=-180/pi*x12, facecolor='white')
    ax.add_artist(p)
    
    p = mpl.patches.Ellipse( (x11+x14/2,x10), x14, x13*x14, angle=180/pi*x12, facecolor='white')
    ax.add_artist(p)

    # draw pupils
    p = mpl.patches.Ellipse( (-x11-x14/2-x15*x14/2, x10), .05, .05, facecolor='black')
    ax.add_artist(p)
    p = mpl.patches.Ellipse( (x11+x14/2-x15*x14/2, x10), .05, .05, facecolor='black')
    ax.add_artist(p)
    
  
#fig = figure(figsize=(11,11))
#for i in range(25):
#    ax = fig.add_subplot(5,5,i+1,aspect='equal')
#    cface(ax, .9, *rand(17))
#    ax.axis([-1.2,1.2,-1.2,1.2])
#    ax.set_xticks([])
#    ax.set_yticks([])

fig.subplots_adjust(hspace=0, wspace=0)
plt.show()
fig = figure(figsize=(10,10))


for i in range(30):
	z=data['z'][i]
	modelMag_u=data['modelMag_u'][i]
	modelMag_r=data['modelMag_r'][i]
	a=np.mean(z/0.2)
	b=np.mean((modelMag_u-20)/3.0)
	c=np.mean((modelMag_u-modelMag_r)/3.0)
	ax = fig.add_subplot(6,6,i+1,aspect='equal')
	cface(ax, a,b,c)
	ax.axis([-1.2,1.2,-1.2,1.2])
	ax.set_xticks([])
	ax.set_yticks([])
	i+=1
