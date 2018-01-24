
% test set file
params.testSetFilename = [params.outFolder 'test_set/testSet.mat'];
% params.testSetFilename = [params.outFolder 'test_set/testSet_diff.mat'];

%% run cubic spline completions

% load test set
data = load(params.testSetFilename,'testSet');
testCurves = data.testSet.curves;


% complete test curves curves
testCompletions = cell(size(testCurves));
for i = 1:numel(testCurves)
%     i
    
    % get inducers
    p1 = testCurves{i}.pts(testCurves{i}.p1, :);
    p2 = testCurves{i}.pts(testCurves{i}.p2, :);
    or1 = testCurves{i}.p1Or;
    or2 = testCurves{i}.p2Or;
    
    % complete curve
    c = completeCurveCS(p1, or1, p2, or2);
    testCompletions{i} = c;
end

% save completions
save([params.outFolder 'test_set/testSet_completions_spline.mat'], 'testCompletions');
% save([params.outFolder 'test_set/testSet_diff_completions_spline.mat'], 'testCompletions');

