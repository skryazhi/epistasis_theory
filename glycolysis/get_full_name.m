function fn = get_full_name(dict, sn)

fn = dict{ strcmp( dict(:,1) , sn ), 2 } ;