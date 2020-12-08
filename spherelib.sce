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

// Converts a decimal to degrees and minutes.
function dm=deg_to_dm(deg)
    dm(1) = floor(deg)
    dm(2) = (deg - dm(1)) * 60
endfunction

// Extension of atan(y/x), which adjusts for quadrants.
// Result is in radians.
function [theta]=atan2(y, x)
    if x >= 0 & y >= 0 then
        theta = atan(abs(y/x))           // Quadrant 1 (+, +)
    elseif x < 0 & y >= 0 then
        theta = atan(abs(x/y)) + %pi/2   // Quadrant 2 (-, +)
    elseif x < 0 & y < 0 then
        theta = atan(abs(y/x)) + %pi     // Quadrant 3 (-, -)
    elseif x >= 0 & y < 0 then
        theta = atan(abs(x/y)) + 3*%pi/2 // Quadrant 4 (+, -)
    end
endfunction

// Extension of atan(y/x), which adjusts for quadrants.
// Result is in degrees.
function [theta]=atan2d(y, x)
    theta = atan2(y, x) * rad_to_deg
endfunction

// Convert a (x, y) point (cartesian) tp polar (r, theta) coordinates.
// Theta is in degrees.
function p=xy_to_polar(xy)
    p(1) = sqrt(xy(1)^2 + xy(2)^2)
    p(2) = atan2d(xy(2), xy(1))
endfunction

// Converts polar (r, theta) coordinates to cartesian (x, y) coordinates.
// Theta is in degrees.
function xy=polar_to_xy(p) 
    xy(1) = p(1) * cosd(p(2))
    xy(2) = p(1) * sind(p(2))
endfunction

// Converts cylindrical (r, theta, z) coordinates to
// spherical (r, theta, phi) coordinates.
// Theta, and phi, are in degrees. 
function xyz=cyl_to_xyz(p)
    xyz(1) = p(1) * cosd(p(2))
    xyz(2) = p(1) * sind(p(2))
    xyz(3) = p(3)
endfunction

// Converts cylindrical (r, theta, z) coordinates to
// spherical (r, theta, phi) coordinates.
// Theta, and phi, are in degrees.
function sph=cyl_to_sph(cyl)
    xyz = cyl_to_xyz(cyl)
    sph = xyz_to_sph(xyz)
endfunction

// Converts cartesian (x, y, z) coordinates to spherical (r, theta, phi)
// coordinates. Theta, and phi, are in degrees.
function sph=xyz_to_sph(xyz)
    sph(1) = sqrt(xyz(1)^2 + xyz(2)^2 + xyz(3)^2)
    sph(2) = atan2d(xyz(2), xyz(1))
    sph(3) = acosd(xyz(3) / sph(1))
endfunction

// Converts spherical coordinates to cartesian coordinates. Spherical
// coordinates are given in degrees!
function p=sph_to_xyz(sph)
    p(1) = sph(1) * sind(sph(3)) * cosd(sph(2))
    p(2) = sph(1) * sind(sph(3)) * sind(sph(2))
    p(3) = sph(1) * cosd(sph(3))
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

// Corrects the angle based on if it should be acute or obtuse (~acute).
function angle=correct_angle(is_acute, theta)
    angle = theta
    if is_acute & angle > 90 | ~is_acute & angle < 90 then
        angle = 180 - angle
    end
endfunction

// Solve a Right Angled Spherical Triangle (RAST) when a, and c is given.
// The result is a vector with [a, b, c, A, B, C], with all angles in degrees.
function u=solve_RST_ac(a, c)
    // Check if b should be acute with cos c = cos a cos b. Listing all 
    // the possible combinations in a table:
    //
    // cos c = cos a * cos b
    // ---------------------
    // + (T) | + (T) | + (T)   doesn't occur
    // - (F) | + (T) | - (F)
    // - (F) | - (F) | + (T)
    // + (T) | - (F) | - (F)
    //
    // From the table we can see that b is acute when:  c & a | ~c & ~a, which
    // is also true when a == c, but only ~c & ~a occurs which is ~(c & a).
    a_is_acute = a < 90
    c_is_acute = c < 90
    b_is_acute = ~(a_is_acute & c_is_acute)
    // Solve for b with cos c = cos a cos b.
    b = acosd(cosd(c) / cosd(a))
    b = correct_angle(b_is_acute, b)
    // Solve for A with sin A = sin a / sin c.
    A = asind(sind(a) / sind(c))
    A = correct_angle(a_is_acute, A) // if a is acute, then A also is.
    // Solve for B with cos A = sin B cos a.
    B = asind(cosd(A) / cosd(a))
    B = correct_angle(b_is_acute, B) // if b is acute, then B also is.
    u = [a, b, c, A, B, 90]
endfunction

// This is a shorter, but less readable, version of solve_RST_ac.
function u=solve_RST_ac_short(a, c)
    b = correct_angle((a > 90) & (c > 90), acosd(cosd(c) / cosd(a)))
    A = correct_angle( a < 90            , asind(sind(a) / sind(c)))
    B = correct_angle((a > 90) & (c > 90), asind(cosd(A) / cosd(a)))
    u = [a, b, c, A, B, 90]
endfunction



