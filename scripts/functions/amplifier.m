function [y] = amplifier(x, type)

type = upper(type);

    switch type
        case 'IDEAL'
            y = x;
        case 'CLIPPING'
            y = min(abs(x), A0) .* exp(1j * angle(x));
        case 'TWT'
            g_A = (chi_A .* abs(x)) ./ (1 + kappa_A .* abs(x).^2);
            g_phi = (chi_phi .* abs(x).^2) ./ (1 + kappa_phi .* abs(x).^2);
            y = g_A .* exp(1i * (angle(x) + g_phi));
        case 'SOLID_STATE'
            y = abs(x) ./ (1 + (abs(x) / A0).^(2 * rho)).^(1 / (2 * rho)) .* exp(1j * angle(x));
        otherwise
            error('Não existe esse parâmetro.');
    end
end