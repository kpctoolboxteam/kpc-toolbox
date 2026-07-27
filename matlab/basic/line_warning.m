function line_warning(caller, msg, varargin)
% LINE_WARNING(CALLER, MSG, ...) - issue a formatted warning tagged with the
% calling function.
warning('%s: %s', caller, sprintf(msg, varargin{:}));
end
