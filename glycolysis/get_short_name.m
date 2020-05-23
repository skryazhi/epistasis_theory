function sn = get_short_name(dict, fn)

sn = dict{ strcmp( dict(:,2) , fn ), 1};