// Conversion constants
rad_to_deg = 180 / %pi

// Calculates the angle and distance of a loxodrome for two given
// points A and B on a sphere. The result is a vector with [angle, dist].
// All values are in degrees.
function r=loxodrome_latlondeg(lat_a_deg, lon_a_deg, lat_b_deg, lon_b_deg)
    delta_phi = lat_b_deg - lat_a_deg
    delta_L = lon_b_deg - lon_a_deg
    delta_VB = rad_to_deg * ( log( tand(45 + lat_b_deg / 2 ) / ..
                                  tand(45 + lat_a_deg / 2 ) ) )
    k = atand(delta_L / delta_VB)
    if delta_L < 0 | k < 0 then
        k = k + 180
    end
    v = delta_phi * 60 / cosd(k)
    r(1) = k
    r(2) = v
endfunction

// Example usage
loxodrome_latlondeg(52 + 11.3/60, 2 + 27/60, 43 + 27.1/60, -(4 + 36.6/60))
