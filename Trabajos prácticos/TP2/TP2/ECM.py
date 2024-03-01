def ECM(X, X_descomp):
	X_flat = X.flatten()
	X_descomp_flat = X_descomp.flatten()
	n = len(X_flat)
	ecm = 0
	i = 0
	while i < n:
		ecm += (X_flat[i] - X_descomp_flat[i])**2
		i += 1
	ecm = ecm / n

	return ecm
