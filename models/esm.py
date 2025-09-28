# models/esm.py
from __future__ import annotations
from typing import Tuple
import torch
import numpy as np

_ESM_SINGLETON: dict | None = None


def _pick_device() -> torch.device:
    if torch.backends.mps.is_available():
        return torch.device("mps")
    if torch.cuda.is_available():
        return torch.device("cuda")
    return torch.device("cpu")


def load_esm_small() -> Tuple[torch.nn.Module, object, object, torch.device, int]:
    """
    Load tiny ESM2 model (t6_8M) to keep PoC fast.
    Returns: (model, alphabet, batch_converter, device, repr_layer)
    """
    global _ESM_SINGLETON
    if _ESM_SINGLETON is not None:
        return (
            _ESM_SINGLETON["model"],
            _ESM_SINGLETON["alphabet"],
            _ESM_SINGLETON["batch_converter"],
            _ESM_SINGLETON["device"],
            _ESM_SINGLETON["repr_layer"],
        )

    from esm import pretrained  # lazy import

    model, alphabet = pretrained.esm2_t6_8M_UR50D()  # ~8M params
    device = _pick_device()
    model = model.to(device)
    model.eval()
    batch_converter = alphabet.get_batch_converter()
    repr_layer = 6  # last layer for t6
    _ESM_SINGLETON = {
        "model": model,
        "alphabet": alphabet,
        "batch_converter": batch_converter,
        "device": device,
        "repr_layer": repr_layer,
    }
    return model, alphabet, batch_converter, device, repr_layer


@torch.no_grad()
def embed_sequence_mean(seq: str) -> np.ndarray:
    """
    Mean-pool ESM token embeddings (exclude BOS/EOS/pad).
    Returns np.float32 vector shape (d,)
    """
    model, alphabet, batch_converter, device, repr_layer = load_esm_small()
    data = [("query", seq)]
    _, _, toks = batch_converter(data)
    toks = toks.to(device)
    out = model(toks, repr_layers=[repr_layer], return_contacts=False)
    reps = out["representations"][repr_layer]  # (B, L, D)
    # Determine true length from sequence (no specials)
    true_len = len(seq)
    # tokens layout: [BOS] + seq + [EOS] (+ padding)
    token_slice = reps[0, 1 : 1 + true_len, :]
    vec = token_slice.mean(dim=0).detach().to("cpu").float().numpy()
    return vec
