Micromorphic (Mindlin form) homogenization pipeline - run:   main_micromorphic          (quick check: main_micromorphic(true))

Model:  W = 1/2 eps'C_hom eps + 1/2 gamma'C_gamma gamma + eps'C_couple gamma + 1/2 l^2 G |grad psi|^2,  gamma = eps - psi  (C_hom fixed)
C_hom:  P.material.ChomSource in parameters.m:  'homogenize' (default, homogenize2D_PBC, cached)  |  'fixed' (uses P.material.ChomFixed)

Reuse:  every stage is cached with its input; unchanged input -> loaded.   Force:  main_micromorphic(false,{'case2'},{'identification'})  or 'all'
        stages: homog, fullIdent, identification, fullTest, micro.   After changing code (not parameters) raise P.out.version.
Output: results/ (results_quick/ for quick runs): homogenization/, case1/{fullscale,identification,logs,figures}, case2/{fullscale,micromorphic,logs,figures,comparison.csv}
        fields: displacement, strain, stress; bin means = before smoothing, Legendre fits = after smoothing; PNGs at the end (P.out.png).
        P.out.rawFields=false stores only bin means + displacements (much smaller files for fine meshes).
