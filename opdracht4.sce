function [v]=vertex(lat_a, lon_a, lat_b, lon_b)
    
    // bA, bB en Lab bepalen
    bA = lat_a
    bB = lat_b
    Lab = abs(lon_b - lon_a)
    if(Lab > 180)
       Lab = 360 - Lab     
    end

    // hoek GrKgrc bepalen
    tanK = sind(Lab) / (cosd(bA)*tand(bB)-sind(bA)*cosd(Lab))
    K = atand(tanK)
    if(K < 0 && Lab > 0) 
        K = K + 180
    elseif(K > 0 && Lab < 0) 
        K = K + 180
    elseif(K < 0 && Lab < 0) 
        K = K + 360
    end        

    // lengte van vertex bepalen
    Lav = atand(1/(tand(K)*sind(bA)))

    LV = lon_a + Lav
    if(LV > 180)
        LV = LV - 360 // West
    end

    // breedte van de vertex bepalen
    tanBV = tand(bA)/cosd(Lav)
    BV = atand(tanBV)
    
    v(1) = BV
    v(2) = LV
endfunction
