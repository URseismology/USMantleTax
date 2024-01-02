function [vp, vs] = find_closest_vp_vs(lat, long, SchultzPelkum)
    distances = sqrt((lat - SchultzPelkum.latitude).^2 + (long - SchultzPelkum.longitude).^2);
    [~, closest_idx] = min(distances);
    vp = SchultzPelkum.vp(closest_idx);
    vs = SchultzPelkum.vsv(closest_idx);
end