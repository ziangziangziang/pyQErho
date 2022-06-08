import numpy as np

def hpsort_eps(ra:np.ndarray, eps=1e-8):
    '''
    to sort the ra array within error eps
    adapted from Modules/sort.f90
    subroutine hpsort_eps(n, ra, ind, eps)
    in Modules/recvec_subs.f90
    hpsort_eps(ngm_g, g2sort_g, igsrt, eps8)
    '''
    n = len(ra)
    ind = np.arange(n)
    if n < 2:
        return ra, ind      # no need to sort
    l = n // 2 # left
    ir = n - 1 # the border at the right side
    while True:
        if l > 0: # so called hiring phase
            # heapify
            l -= 1
            rra = ra[l]     # root
            iind = ind[l]   # index of root
        else:
            rra = ra[ir]
            iind = ind[ir]
            ra[ir] = ra[0]  # retire the top of the heap into array
            ind[ir] = ind[0]# retire the top of the heap into the index
            ir -= 1
            if ir == 0:
                ra[0] = rra
                ind[0] = iind
                break

        # either in hiring or promotion phase
        i = l
        j = l+l+1
        while j<=ir:
            if j<ir:
                if abs(ra[j] - ra[j+1]) >= eps:
                    if ra[j] < ra[j+1]:
                        j += 1 # j is the larger one
                elif ind[j] < ind[j+1]:
                    j+=1
            # demote rra
            if abs(rra - ra[j]) >= eps:
                if rra<ra[j]:
                    ra[i] = ra[j]
                    ind[i] = ind[j]
                    i = j
                    j = i+i+1
                else:
                    j = ir+1
            else:
                # rra == ra[j] within tolerance
                # then sort by index value
                if iind < ind[j]:
                    ra[i] = ra[j]
                    ind[i] = ind[j]
                    i = j
                    j = i+i+1
                else:
                    # set j tp ternimate loop
                    j = ir+1
        ra[i] = rra
        ind[i] = iind
    return ra, ind

if __name__ == '__main__':
    import time
    import matplotlib.pyplot as plt
    ra = np.random.random(110000)
    start = time.time()
    rb, ind = hpsort_eps(ra)
    print(time.time()-start)
    rc = np.sort(ra, kind='heapsort')
    start = time.time()
    print(time.time()-start)
    plt.plot(rb)
    plt.plot(rc)
    plt.show()


