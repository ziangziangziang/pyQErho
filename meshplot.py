import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np

def explode(data):
    size = np.array(data.shape)*2
    if data.ndim == 3:
        data_e = np.zeros(size-1, dtype=data.dtype)
        data_e[::2, ::2, ::2] = data
    elif data.ndim == 4:
        size[-1] = size[-1]//2
        data_e = np.zeros(size-np.array((1,1,1,0)), dtype=data.dtype)
        data_e[::2, ::2, ::2, ] = data

    return data_e

def linear_alpha(data):
    data_max = data.max()
    data_min = data.min()
    f = lambda x: ((x-data_min)/(data_max-data_min))**2
    return f(data)

def colorfunction(data):
    norm = matplotlib.colors.Normalize(vmin=data.min(), vmax=data.max(), clip=True)
    mapper = cm.ScalarMappable(norm=norm, cmap=cm.Plasma)
    rgba_array = np.zeros(data.shape, dtype=np.float)
    x, y, z = data.shape
    for i in range(x):
        for j in range(y):
            for k in range(z):
                rgba_array[i, j, k] = mapper.to_rgba(data[i, j, k])
    return rgba_array

def mesh_plot(data, output_name=None, filtering='none', title='none', a=np.identity(3)):
    # uncomment if to plot 3D
    data_e = explode(data)
    print(data.shape)

    #fcolors = colorfunction(data)
    fcolors = plt.cm.plasma(data)
    alpha_dens = linear_alpha(data)
    fcolors[:, :, :, 3] = alpha_dens

    #print(fcolors[0,0,0]*255)
    #print(fcolors[5,5,5]*255)
    norm = matplotlib.colors.Normalize(vmin=0, vmax=data.max())
    #filled = np.ones(data.shape)
    if filtering == 'average':
        filled = np.where(data>=data.mean(), 1, 0)
    elif filtering == 'reverse average':
        filled = np.where(data<=data.mean(), 1, 0)
    elif filtering == 'none':
        filled = np.ones_like(data)
    elif filtering == 'half':
        filled = np.zeros_like(data)
        num_x, num_y, num_z =  filled.shape
        filled[:, num_y//2:, :] = 1
    x, y, z = np.indices(np.array(data_e.shape)+1).astype(float)
    x = x/(data_e.shape[0]+1)
    y = y/(data_e.shape[1]+1)
    z = z/(data_e.shape[2]+1)
    x[0::2, :, :] += 0.05/(data_e.shape[0]+1)
    y[:, 0::2, :] += 0.05/(data_e.shape[1]+1)
    z[:, :, 0::2] += 0.05/(data_e.shape[2]+1)
    x[1::2, :, :] += 0.95/(data_e.shape[0]+1)
    y[:, 1::2, :] += 0.95/(data_e.shape[1]+1)
    z[:, :, 1::2] += 0.95/(data_e.shape[2]+1)

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    #ax.set_box_aspect(data.shape)
    #ax.set_aspect('auto')
    #print(data.mean())
    #print(data.min())
    #print(data.max())
    ax.set_xlabel('x $\AA$')
    ax.set_ylabel('y $\AA$')
    ax.set_zlabel('z $\AA$')
    #ax.set_aspect('auto')
    ax.set_title(title)
    print('calculation done, rendering voxels....\n This may take lots of time')
    a = np.array(a)
    ax.quiver(0, 0, 0, a[0][0],  a[0][1],  a[0][2], arrow_length_ratio = 0.1,color='k', alpha=0.6)
    ax.quiver(0, 0, 0, a[1][0],  a[1][1],  a[1][2], arrow_length_ratio = 0.1,color='k', alpha=0.6)
    ax.quiver(0, 0, 0, a[2][0],  a[2][1],  a[2][2], arrow_length_ratio = 0.1,color='k', alpha=0.6)
    #ax.plot(0, 0, 0, linestyle=None, marker='o')
    #ax.quiver(0, 0, 0, 0,  0,  1, arrow_length_ratio = 0.1,color='k', alpha=0.6)
    #ax.quiver(0, 0, 0, 0,  1,  0, arrow_length_ratio = 0.1,color='k', alpha=0.6)
    #ax.quiver(0, 0, 0, 1,  0,  0, arrow_length_ratio = 0.1,color='k', alpha=0.6)

    x,y,z = np.matmul(np.array([x,y,z]).transpose((1,2,3,0)), a).transpose(3,0,1,2)
    max_range = np.array([x.max()-x.min(), y.max()-y.min(), z.max()-z.min()])

    ax.set_box_aspect(max_range)
    ax.voxels(x, y, z, explode(filled),\
            facecolors=explode(fcolors))
    m = cm.ScalarMappable(cmap = plt.cm.plasma, norm=norm)
    m.set_array([])
    plt.colorbar(m, label='charge density', ticks=[data.max(), data.min()])
    if output_name:
        plt.savefig(output_name, dpi=300)
    else:
        plt.show()
    plt.cla()
    plt.close(fig)
    

