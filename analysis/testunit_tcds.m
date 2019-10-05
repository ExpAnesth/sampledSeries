function tests=testunit_tcds
% tests for function tcds - only very first, simple tests implemented
tests=functiontests(localfunctions);
end

function setupOnce(testCase)
% test input #1
testCase.TestData.d1=[zeros(10,1); ones(3,1); zeros(10,1)];
% test input #2
testCase.TestData.d2=[zeros(10,1); ones(1,1); zeros(5,1); ones(1,1)];
figure(1), clf
plot(testCase.TestData.d1,'o-');
hold on 
plot(testCase.TestData.d2,'o-')
niceyax
end

function test01(testCase)
% too high threshold should result in empty tsl
tsl = tcds(testCase.TestData.d1, 'idx', 2);
assert(isempty(tsl), 'nonempty output with too high threshold')
end

function test02(testCase)
tsl = tcds(testCase.TestData.d1, 'idx', 0.2);
assert(tsl==11, 'bad tsl')
end

