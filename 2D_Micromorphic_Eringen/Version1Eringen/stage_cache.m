function [D, computed] = stage_cache(file, key, reuse, fcn, logFile)
% =========================================================================
%  STAGE CACHE      [D, computed] = stage_cache(file, key, reuse, fcn, logFile)
% =========================================================================
%  Returns the data D stored in 'file' if the file was made with exactly the same input 'key' (and reuse is true).
%  Otherwise runs D = fcn(), saves D and key to 'file' (MAT-file v7.3) and returns D.
%      key      : struct with EVERYTHING the result depends on (compared with isequal)
%      reuse    : false forces a new calculation
%      logFile  : (optional) the screen output of fcn is written to this text file (only when fcn is run)
%      computed : true if fcn was run, false if D was loaded
%  Only the (small) key is read to decide; the (large) data are loaded only when the key matches.
% =========================================================================
if nargin < 5, logFile = ''; end
computed = false;
if reuse && exist(file, 'file') == 2
    try
        s = load(file, 'key');
        if isequal(s.key, key)
            t = load(file, 'D');  D = t.D;
            fprintf('  [cache] loaded  %s  (input unchanged)\n', file);
            return
        end
        fprintf('  [cache] %s was made with other input -> recomputing\n', file);
    catch
        fprintf('  [cache] %s is not readable -> recomputing\n', file);
    end
end
lg = [];
if ~isempty(logFile), lg = start_log(logFile); end %#ok<NASGU>
D = fcn();
computed = true;
folder = fileparts(file);
if ~isempty(folder) && exist(folder, 'dir') ~= 7, mkdir(folder); end
save(file, 'D', 'key', '-v7.3');
fprintf('  [cache] saved   %s\n', file);
clear lg                                                                    % closes the log
end
