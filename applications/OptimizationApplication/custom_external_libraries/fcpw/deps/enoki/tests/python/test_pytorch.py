import enoki as ek
import numpy as np
import pytest
import torch

class EnokiAtan2(torch.autograd.Function):
    """PyTorch function example from the documentation."""
    @staticmethod
    def forward(ctx, arg1, arg2):
        ctx.in1 = ek.FloatD(arg1)
        ctx.in2 = ek.FloatD(arg2)
        ek.set_requires_gradient(ctx.in1, arg1.requires_grad)
        ek.set_requires_gradient(ctx.in2, arg2.requires_grad)
        ctx.out = ek.atan2(ctx.in1, ctx.in2)
        out_torch = ctx.out.torch()
        ek.cuda_malloc_trim()
        return out_torch

    @staticmethod
    def backward(ctx, grad_out):
        ek.set_gradient(ctx.out, ek.FloatC(grad_out))
        ek.FloatD.backward()
        result = (ek.gradient(ctx.in1).torch()
                  if ek.requires_gradient(ctx.in1) else None,
                  ek.gradient(ctx.in2).torch()
                  if ek.requires_gradient(ctx.in2) else None)
        del ctx.out, ctx.in1, ctx.in2
        ek.cuda_malloc_trim()
        return result


def test01_set_gradient():
    a = ek.FloatD(42, 10)
    ek.set_requires_gradient(a)

    with pytest.raises(TypeError):
        grad = ek.FloatD(-1, 10)
        ek.set_gradient(a, grad)

    grad = ek.FloatC(-1, 10)
    ek.set_gradient(a, grad)
    assert np.allclose(grad.numpy(), ek.gradient(a).numpy())

    # Note: if `backward` is not called here, test03 segfaults later.
    # TODO: we should not need this, there's most likely some missing cleanup when `a` is destructed
    ek.FloatD.backward()
    del a, grad


def test02_array_to_torch():
    a = ek.FloatD(42, 10)
    a_torch = a.torch()
    assert isinstance(a_torch, torch.Tensor)
    a_torch += 8
    a_np = a_torch.cpu().numpy()
    assert isinstance(a_np, np.ndarray)
    assert np.allclose(a_np, 50)


def test03_pytorch_function():
    enoki_atan2 = EnokiAtan2.apply

    y = torch.tensor(1.0, device='cuda')
    x = torch.tensor(2.0, device='cuda')
    y.requires_grad_()
    x.requires_grad_()

    o = enoki_atan2(y, x)
    o.backward()
    assert np.allclose(y.grad.cpu(), 0.4)
    assert np.allclose(x.grad.cpu(), -0.2)
