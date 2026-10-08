function c = start_log(file)
% Starts a diary (screen log) in 'file', replacing an older one.  The diary is closed when the returned object is cleared
% (clear c) or when the function that called start_log ends (also after an error).
folder = fileparts(file);
if ~isempty(folder) && exist(folder, 'dir') ~= 7, mkdir(folder); end
if exist(file, 'file') == 2, delete(file); end
diary('off');
diary(file);
c = onCleanup(@() diary('off'));
end
