import numpy as np

def get_spectral_vector(t, X):
    dt = t[1] - t[0] # in ms
    N = X.shape[0]
    T = N * dt  # in
    df = 1 / T
    fNQ = 1 / dt / 2
    f = np.arange(0, fNQ,df)
    f = f[1:]

    X_f = np.zeros([len(f), X.shape[1]])
    for k in range(X.shape[1]):
        x_f = np.fft.fft(X[:, k] - np.mean(X[:, k]))
        S_f = 2 * dt ** 2 / T * (x_f * x_f.conj())
        S_f = np.real(S_f[1:int(N / 2)])
        X_f[:, k] = S_f

    return f, X_f