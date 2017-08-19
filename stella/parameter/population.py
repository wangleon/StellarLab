import math
import numpy as np

def get_pop_prob(U, V, W, norm=False):
    '''Get relative probability of thin disk, thick disk and halo based on
    kinetic parameters

    Parameters
    ----------
    U : *float*
        Galactic velocity component *U* in km/s relative to the LSR
    V : *float*
        Galactic velocity component *V* in km/s relative to the LSR
    W : *float*
        Galactic velocity component *W* in km/s relative to the LSR
    norm : *bool*
        Normalized the probabilites to 1 if `True`

    Returns
    --------
    prob_lst : *tuple*
        Relative probabilities of thin disk, thick disk and halo

    Notes
    ------
    Calculate the relative probability of stellar population (thin disk/thick
    disk/halo) using the Bayesian method proposed by Bensby et al. 2003. These
    authors assumed that the Galactic space velocities (*U*:sub:`LSR`,
    *V*:sub:`LSR` and *W*:sub:`LSR`) of different popolutaions have Gaussian
    distributions with different characteristic velocity dispersions.
    The probabilities can be calculated as:

    .. math::

        f(U, V, W) = k\exp\left(
                    -\\frac{U_\mathrm{LSR}^2}{2\sigma_U^2}
                    -\\frac{(V_\mathrm{LSR}-V_\mathrm{asym})^2}{2\sigma_V^2}
                    -\\frac{W_\mathrm{LSR}^2}{2\sigma_W^2}
                    \\right)

    where *V*:sub:`asym` is the asymmetric velocity drift, and

    .. math::

        k = \\frac{1}{(2\pi)^{3/2}\sigma_U\sigma_V\sigma_W}

    is the normalize factor.
    
    For a given star, determination of its beloging population should take
    account of the local number densities of different populations.
    In the solar neighbourhood, about 94% of stars belong to the thin disk, and
    6% belong to the thick disk. Therefore, the probability of belonging to a
    specific population (c) is

    .. math::

        p_\mathrm{c} = X_\mathrm{c}f_\mathrm{c}(U, V, W)

    The relative probabilities for the thick-disk-to-thin-disk (TD/D) and
    thick-disk-to-halo (TD/H) membership are

    .. math::

        \mathrm{TD}/\mathrm{D} = \\frac{X_\mathrm{TD}}{X_\mathrm{D}} \\frac{f_\mathrm{TD}}{f_\mathrm{D}},
        \quad\quad
        \mathrm{TD}/\mathrm{H} = \\frac{X_\mathrm{TD}}{X_\mathrm{H}} \\frac{f_\mathrm{TD}}{f_\mathrm{H}}

    The number densities in solar neighbourhood, velocity dispersions and
    asymmetric velocity drifts for thin disk, thick disk and halo stars are
    listed as below.

    .. |kms| replace:: km s\ :sup:`−1`

    .. csv-table:: Characteristic velocity dispersions (Table 1 of Bensby et al. 2003)
        :header: "Population", "*X*", "*σ*:sub:`U` (|kms|)", "*σ*:sub:`V` (|kms|)", "*σ*:sub:`W` (|kms|)", "*V*:sub:`asym` (|kms|)"
        :widths: 30, 20, 25, 25, 25, 30

        "Thin disk (D)",   "0.94",   "35",  "20", "16", "−15"
        "Thick disk (TD)", "0.06",   "67",  "38", "35", "−46"
        "Halo (H)",        "0.0015", "160", "90", "90", "−220"

    Examples
    --------
    Calculate the probability of belonging population of `Arcturus
    <http://simbad.u-strasbg.fr/simbad/sim-id?Ident=Arcturus>`_ (α Boo,
    HIP 69673).
    The heliocentric radial velocity is −5.19 |kms| (SIMBAD). The UVW velocities
    of sun relative to LSR are (−9.6, 14.6, 9.3) |kms| according to `Reid et al.
    2014 <http://adsabs.harvard.edu/abs/2014ApJ...783..130R>`_.

    .. code-block:: python
        
        from stella.catalog.find_catalog import find_HIP2
        from stella.kinetics.orbit import compute_UVW
        from stella.parameter.population import get_pop_prob
        
        item = find_HIP2(69673)
        rv = -5.19                 # from SIMBAD
        ra, dec = item['RAdeg'], item['DEdeg']
        parallax = item['Plx']
        pm = (item['pmRA'], item['pmDE'])
        u, v, w = compute_UVW(ra=ra,dec=dec,parallax=parallax,rv=rv,pm=pm)
        solar_uvw = (-9.6, 14.6, 9.3) # from Reid et al. 2014, ApJ, 783, 130
        u_lsr = u + solar_uvw[0]
        v_lsr = v + solar_uvw[1]
        w_lsr = w + solar_uvw[2]
        probs = get_pop_prob(u_lsr, v_lsr, w_lsr, norm=True)
        # probs = (0.015746971553225683, 0.98174375842748596, 0.002509270019288334)

    References
    -----------
    * `Bensby et al., 2003, A&A, 410, 527 <http://adsabs.harvard.edu/abs/2003A&A...410..527B>`_
    * `Reid et al. 2014, ApJ, 783, 130 <http://adsabs.harvard.edu/abs/2014ApJ...783..130R>`_

    '''

    X = {'thin':0.94, 'thick':0.06, 'halo':0.0015}
    e_U = {'thin': 35.0, 'thick': 67.0, 'halo': 160.0}
    e_V = {'thin': 20.0, 'thick': 38.0, 'halo':  90.0}
    e_W = {'thin': 16.0, 'thick': 35.0, 'halo':  90.0}
    Vasym = {'thin': -15.0, 'thick': -46.0, 'halo': -220.0}

    prob_lst = []
    for pop in ['thin', 'thick', 'halo']:
        k = 1.0/(math.pow(2*math.pi,1.5)*e_U[pop]*e_V[pop]*e_W[pop])
        prob = k*math.exp( - U**2/2.0/e_U[pop]**2
                           - (V-Vasym[pop])**2/2.0/e_V[pop]**2
                           - W**2/2.0/e_W[pop]**2 )
        prob_lst.append(prob*X[pop])

    if norm:
        prob_lst = np.array(prob_lst)
        prob_lst /= prob_lst.sum(dtype=np.float64)

    return tuple(prob_lst)
