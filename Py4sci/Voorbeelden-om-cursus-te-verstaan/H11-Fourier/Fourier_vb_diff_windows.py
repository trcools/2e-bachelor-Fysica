import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
from scipy.fft import fft, fftshift

def compare_windows(M=51, nfft=2047):
    win_fns = {
        "cosine": signal.windows.cosine,
        "bartlett": signal.windows.bartlett,
        "blackman": signal.windows.blackman,
    }

    plt.close("windowing_compare")
    fig, (axw, axf) = plt.subplots(1, 2, figsize=(10, 4), num="windowing_compare")

    for name, fn in win_fns.items():
        w = fn(M)
        axw.plot(w, label=name)

        A = fft(w, nfft) / (len(w) / 2.0)
        freq = np.linspace(-0.5, 0.5, len(A))
        resp = 20 * np.log10(np.abs(fftshift(A / np.abs(A).max())))
        axf.plot(freq, resp, label=name)

    axw.set_title("Windows in tijdsdomein")
    axw.set_xlabel("sample")
    axw.set_ylabel("amplitude")
    axw.legend()

    axf.set_title("Frequency response (genormaliseerd)")
    axf.set_xlabel("genormaliseerde frequentie [cycles/sample]")
    axf.set_ylabel("magnitude [dB]")
    axf.set_ylim(-120, 5)
    axf.legend()

    plt.tight_layout()

compare_windows()
plt.show()