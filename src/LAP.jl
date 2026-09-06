# Prepare the compatible source and its k-derivative
# according to the selected source extension.
function prepare_source_extension(fem, source_loads, Φres, source_mode; atol=1e-6, rtol=1e-4)
    b = source_loads.b
    b_der = source_loads.b_der

    if source_mode === :original

        check = check_source_orthogonality(source_loads, Φres; atol=atol, rtol=rtol)

        if !check.is_orthogonal
            throw(
                ArgumentError(
                    "source_mode = :original requires the original source " *
                    "to be orthogonal to the resonant eigenspace. " *
                    "Use :frozen_correction or :projected_correction instead.",
                ),
            )
        end

        return (;
            b_effective=b,
            b_der_effective=b_der,
            source_mode,
            source_projection_coefficients=nothing,
            derivative_projection_coefficients=nothing,
        )

    elseif source_mode === :frozen_correction

        projection_b = project_load_to_range(fem, b, Φres)

        # The correction determined at k* is frozen with respect to k.
        # Hence its k-derivative is zero.
        b_effective = projection_b.projected
        b_der_effective = b_der

        return (;
            b_effective,
            b_der_effective,
            source_mode,
            source_projection_coefficients=projection_b.coefficients,
            derivative_projection_coefficients=nothing,
        )

    elseif source_mode === :projected_correction

        projection_b = project_load_to_range(fem, b, Φres)

        projection_db = project_load_to_range(fem, b_der, Φres)

        b_effective = projection_b.projected
        b_der_effective = projection_db.projected

        return (;
            b_effective,
            b_der_effective,
            source_mode,
            source_projection_coefficients=projection_b.coefficients,
            derivative_projection_coefficients=projection_db.coefficients,
        )

    else
        throw(
            ArgumentError(
                "Unknown source_mode = $source_mode. " *
                "Use :original, :frozen_correction, or :projected_correction.",
            ),
        )
    end
end


# Compute a particular solution of the singular Helmholtz system.
function compute_particular_solution(
    fem,
    source_loads,
    k,
    Φres;
    source_mode,
    restart=50,
    maxiter=25000,
    reltol=1e-10,
    atol=1e-6,
    rtol=1e-4,
)
    source_data =
        prepare_source_extension(fem, source_loads, Φres, source_mode; atol=atol, rtol=rtol)

    b_effective = source_data.b_effective
    b_der_effective = source_data.b_der_effective

    compatibility_residual = norm(Φres' * b_effective)

    compatibility_scale = max(norm(Φres) * norm(b_effective), 1.0)

    if compatibility_residual > atol + rtol * compatibility_scale

        error(
            "The effective source is not sufficiently orthogonal " * "to the resonant eigenspace.",
        )
    end

    A = fem.K - k^2 * fem.M

    T = promote_type(eltype(A), eltype(b_effective))

    u_special = zeros(T, length(b_effective))

    u_special, history =
        gmres!(u_special, A, b_effective; restart=restart, maxiter=maxiter, reltol=reltol, log=true)

    linear_residual = norm(A * u_special - b_effective) / max(norm(b_effective), eps(Float64))

    return (;
        u_special,
        history,
        b_effective,
        b_der_effective,
        source_mode,
        source_projection_coefficients=source_data.source_projection_coefficients,
        derivative_projection_coefficients=source_data.derivative_projection_coefficients,
        compatibility_residual,
        linear_residual,
    )
end


# Compute the LAP correction and the final LAP solution.
function compute_lap_solution(fem, particular, Φres, k)
    u_special = particular.u_special
    b_der_effective = particular.b_der_effective

    G = Φres' * fem.M * Φres

    derivative_coupling = Φres' * b_der_effective

    rhs = -derivative_coupling - 2k * (Φres' * fem.M * u_special)

    coefficients = (2k * G) \ rhs

    u_lap = u_special + Φres * coefficients

    constraint_residual = norm(derivative_coupling + 2k * (Φres' * fem.M * u_lap))

    return (; u_lap, coefficients, derivative_coupling, constraint_residual)
end


# Compute the particular solution and the final LAP solution.
function solve_lap(
    fem,
    source_loads,
    Φres,
    k;
    source_mode,
    restart=50,
    maxiter=25000,
    reltol=1e-10,
    atol=1e-6,
    rtol=1e-4,
)
    particular = compute_particular_solution(
        fem,
        source_loads,
        k,
        Φres;
        source_mode=source_mode,
        restart=restart,
        maxiter=maxiter,
        reltol=reltol,
        atol=atol,
        rtol=rtol,
    )

    lap = compute_lap_solution(fem, particular, Φres, k)

    return (; particular, lap)
end
