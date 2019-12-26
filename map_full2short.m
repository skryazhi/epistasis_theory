function snl = map_full2short(dict, fnl)

snl = cellfun( @(fn) get_short_name(dict,fn) , fnl, 'UniformOutput', false);
