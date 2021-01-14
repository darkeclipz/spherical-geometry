// Calculates the great circle (GRC) distance (in units of the radius) between 
// two points on a sphere, given with latitude/longitude, and a radius. 
// The cosine rule is used to find the distance.
function [distance]=GRC_distance_latlondeg(lat_a_deg, lon_a_deg, ..
                                           lat_b_deg, lon_b_deg, ..
                                           radius)
    // cos v = sin(phi_a)*sin(phi_b) + cos(phi_a)*cos(phi_b)*cos(delta_L)
    delta_L = abs(lon_b_deg - lon_a_deg)
    cosv = sind(lat_a_deg)*sind(lat_b_deg) + ..
           cosd(lat_a_deg)*cosd(lat_b_deg)*cosd(delta_L)
    v = acosd(cosv)
    distance = v * 60
endfunction

// Example usage
GRC_distance_latlondeg(13 + 18/60, 124 + 47/60, 21 + 25/60, -(157 + 26/60))
