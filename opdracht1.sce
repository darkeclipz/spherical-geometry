// Conversion constants
deg_to_rad = %pi / 180
rad_to_deg = 180 / %pi
deg_to_nm = 60
nm_to_km = 1.852
deg_to_km = deg_to_nm * nm_to_km

// Converts days, minutes, seconds to degrees.
// Note that d S (south) is entered as -d, and for
// d W (west), it is also entered as -d.
// Example: 5, 2' 11'' W = dms_to_deg(-5, 2, 11)
// A positive result indicates N, E, and a negative
// result indicates S, W.
function deg=dms_to_deg(d, m, s)
    deg = sign(d) * (abs(d) + abs(m) / 60 + abs(s) / (60 * 60))
endfunction

// Converts a latitude and longitude (in degrees) into 
// spherical coordinates (theta, phi) in degrees.
function sph=latlon_to_sph(lat_deg, lon_deg, radius)
    sph(1) = radius
    // Calculate theta
    if lon_deg >= 0 then
        sph(2) = lon_deg
    else
        sph(2) = 360 - abs(lon_deg)
    end
    // Calculate phi
    if lat_deg >= 0 then
        sph(3) = 90 - lat_deg
    else
        sph(3) = 180 - abs(lat_deg)
    end
endfunction

// Converts spherical coordinates to cartesian coordinates. Spherical
// coordinates are given in degrees!
function p=sph_to_xyz(sph)
    p(1) = sph(1) * sind(sph(3)) * cosd(sph(2))
    p(2) = sph(1) * sind(sph(3)) * sind(sph(2))
    p(3) = sph(1) * cosd(sph(3))
endfunction

// Calculates the angle between two vectors in degrees.
function [angle_deg]=vangle(u, v)
    angle_deg = acosd(u'*v / (norm(u) * norm(v)))
endfunction

// Calculates the great circle (GRC) distance (in units of the radius) between 
// two points on a sphere, given with latitude/longitude, and a radius. 
// The dot-product formula is used to calculate the distance in units
// of the radius.
function [distance]=GRC_distance_latlondeg(lat_a_deg, lon_a_deg, ..
                                           lat_b_deg, lon_b_deg, ..
                                           radius)
    sph_a = latlon_to_sph(lat_a_deg, lon_a_deg, radius)
    sph_b = latlon_to_sph(lat_b_deg, lon_b_deg, radius)
    xyz_a = sph_to_xyz(sph_a)
    xyz_b = sph_to_xyz(sph_b)
    angle_deg = vangle(xyz_a, xyz_b)
    distance = angle_deg * deg_to_rad * radius
endfunction
