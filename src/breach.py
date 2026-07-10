"""Dam-breach outflow hydrographs (prescribed, empirical — not structural).

Two breach modes:
  * "froehlich": Froehlich (2008) regression for average breach width and
    formation time; trapezoidal breach grows linearly (invert drops, width
    opens) over the formation time.
  * "sudden": a fixed-width rectangular gap (e.g. 1/3 of crest length) opens
    to full depth instantaneously — the "masonry monolith slides out" case.

Outflow is a broad-crested weir over the instantaneous breach geometry,
coupled to a level-pool reservoir with a V-shaped (quadratic) stage-storage
curve calibrated so the initial pool matches the scenario volume. All
elevations here are heads above the breach bottom (dam base = 0).
"""
from dataclasses import dataclass, field

import numpy as np

G = 9.81
CW_RECT = 1.7   # SI broad-crested weir coefficient, rectangular part [m^0.5/s]
CW_TRI = 1.35   # SI coefficient for the triangular side areas


@dataclass
class BreachSpec:
    name: str
    volume_m3: float        # reservoir volume at failure
    head_m: float           # initial water surface above breach bottom
    crest_len_m: float
    mode: str = "froehlich"       # or "sudden"
    k0: float = 1.0               # Froehlich failure-mode factor (1.3 overtopping)
    side_slope: float = 1.0       # trapezoid side slope m (H:V), froehlich mode
    sudden_frac: float = 1.0 / 3.0  # fraction of crest removed in sudden mode
    # derived
    b_avg: float = field(init=False, default=0.0)
    t_form: float = field(init=False, default=0.0)

    def __post_init__(self):
        if self.mode == "froehlich":
            self.b_avg = 0.27 * self.k0 * self.volume_m3**0.32 * self.head_m**0.04
            self.b_avg = min(self.b_avg, self.crest_len_m)
            self.t_form = 63.2 * np.sqrt(self.volume_m3 / (G * self.head_m**2))
        else:  # sudden
            self.b_avg = self.crest_len_m * self.sudden_frac
            self.t_form = 0.0


def hydrograph(spec: BreachSpec, duration_s: float = 86400.0, dt: float = 1.0,
               out_dt: float = 60.0):
    """Integrate the level-pool / weir ODE.

    Returns (t_out, q_out) sampled every `out_dt` seconds, plus peak info.
    """
    # V-shape stage-storage: V(eta) = 0.5 * a * eta^2, calibrated to the pool.
    a_coef = 2.0 * spec.volume_m3 / spec.head_m**2

    v = spec.volume_m3
    n_out = int(duration_s / out_dt) + 1
    t_out = np.arange(n_out) * out_dt
    q_out = np.zeros(n_out)

    q_peak, t_peak = 0.0, 0.0
    t = 0.0
    k = 0
    while t <= duration_s:
        eta = np.sqrt(max(2.0 * v / a_coef, 0.0))  # pool surface above base
        if spec.mode == "froehlich":
            f = min(t / spec.t_form, 1.0) if spec.t_form > 0 else 1.0
            z_inv = spec.head_m * (1.0 - f)          # invert drops crest->base
            b_bot_final = max(spec.b_avg - spec.side_slope * spec.head_m, 0.0)
            b_bot = b_bot_final * f
            m = spec.side_slope
        else:
            z_inv, b_bot, m = 0.0, spec.b_avg, 0.0

        head = eta - z_inv
        if head > 0.0:
            q = CW_RECT * b_bot * head**1.5 + CW_TRI * m * head**2.5
            q = min(q, v / dt)  # cannot release more than what remains
        else:
            q = 0.0

        if q > q_peak:
            q_peak, t_peak = q, t

        while k < n_out and t_out[k] <= t:
            q_out[k] = q
            k += 1

        v = max(v - q * dt, 0.0)
        t += dt

    return t_out, q_out, q_peak, t_peak


def mullaperiyar_spec(pool_ft: int, mode: str) -> BreachSpec:
    """The two Mullaperiyar pools used by the scenarios."""
    if pool_ft == 142:
        vol, head = 380e6, 43.3
    elif pool_ft == 152:
        vol, head = 443e6, 46.3
    else:
        raise ValueError(pool_ft)
    spec = BreachSpec(
        name=f"mullaperiyar_{pool_ft}ft_{mode}",
        volume_m3=vol, head_m=head, crest_len_m=365.7, mode=mode, k0=1.0,
    )
    # experiment knob: override the Froehlich formation time (e.g. the
    # IIT Roorkee study assumed a 12-minute breach)
    import os
    tf_min = os.environ.get("MULLA_BREACH_TF_MIN")
    if tf_min is not None and spec.mode == "froehlich":
        spec.t_form = float(tf_min) * 60.0
        spec.name += f"_tf{tf_min}min"
    return spec


def cheruthoni_spec() -> BreachSpec:
    """Cascade failure of Cheruthoni (proxy for the Idukki cluster).

    Treated as overtopping-induced (K0=1.3) when the Mullaperiyar surge
    arrives. Live storage dwarfs Mullaperiyar — that is the point.
    """
    return BreachSpec(
        name="cheruthoni_cascade",
        volume_m3=1460e6, head_m=138.0, crest_len_m=650.0,
        mode="froehlich", k0=1.3,
    )


if __name__ == "__main__":
    for spec in [
        mullaperiyar_spec(142, "froehlich"),
        mullaperiyar_spec(152, "sudden"),
        cheruthoni_spec(),
    ]:
        t, q, qp, tp = hydrograph(spec)
        vol_out = np.trapezoid(q, t)
        print(
            f"{spec.name}: B_avg={spec.b_avg:.0f} m  t_form={spec.t_form/3600:.2f} h  "
            f"Qpeak={qp:,.0f} m3/s at t={tp/60:.1f} min  "
            f"released={vol_out/1e6:.0f} Mm3 / {spec.volume_m3/1e6:.0f} Mm3"
        )
