"""Physics as a signal: residuals, PINN machinery and sensitivities.

Three distinct notions of "residual" live here, and confusing them is the usual
mistake:

``solver_residuals``
    ``ResidualEvaluator`` assembles the real PDE residual of a predicted field
    through the solver's own builder. Cheap, **not** differentiable - a *score*
    for query strategies, callbacks and validation.
``physics_informed``
    SymPy strong-form residuals evaluated by ``physicsnemo.sym``'s
    ``PhysicsInformer`` as differentiable **training loss terms**. Approximate:
    it is the PDE, not the discretization.
``differentiable_residual``
    The exact discrete residual through the real FEM assembly, wrapped as a
    ``torch.autograd.Function`` (forward = ``BuildRHS``, backward = the
    consistent tangent's transpose). Exact and differentiable, and the most
    expensive.

``sensitivity_utils`` sits alongside them: cheap surrogate ``dJ/dx`` by autograd,
exact adjoint parameter sensitivities, and the discretely exact shape gradient
at every node from one pass over the mesh.

``diffusion_residual_operator`` puts the third notion to a different use: the
same assembly, wrapped as the observation operator of a diffusion posterior
sampling guidance term, so the physics steers a generative model at SAMPLING
time rather than grading it during training.
"""
