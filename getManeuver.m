function maneuver_nb = getManeuver(maneuver)
    switch maneuver
        case 'none'
            maneuver_nb = 1;
        case 'DLC'
            maneuver_nb = 2;
        case 'split-mu'
            maneuver_nb = 3;
        otherwise
            error('Check the value of varaible maneuver')
    end
end

