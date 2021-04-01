// Copyright 2021 Frolova Olga
#include <math.h>
#include <omp.h>
#include <random>
#include <vector>
#include <algorithm>
#include <iostream>
#include <string>
#include "../../../modules/task_1/frolova_o_radix_sort_batcher_merge/radix_sort_batcher_merge.h"

using namespace std;

std::vector<std::pair<int, int>> comps;

std::vector<double> getRandomVector(int size) {
    std::vector<double> vec(size);
    std::mt19937 gen;
    gen.seed(static_cast<unsigned char>(time(0)));
    for (int i = 0; i < size; i++) {
        vec[i] = gen() / 10000;
    }
    return vec;
}

int countRadix(double num) {
    int res = 0;
    while (static_cast<int>(num) > 0) {
        num = num / 10;
        res++;
    }
    return res;
}

std::vector<double> checkVector(std::vector<double> vec) {
    int j = 0;
    double tmp = 0;
    int lengh = static_cast<int>(vec.size());
    for (int i = 0; i < lengh; i++) {
        j = i;
        for (int k = i; k < lengh; k++) {
            if (vec[j] > vec[k]) {
                j = k;
            }
        }
        tmp = vec[i];
        vec[i] = vec[j];
        vec[j] = tmp;
    }
    return vec;
}

double maxVal(std::vector<double> vec) {
    int size = static_cast<int>(vec.size());
    double res = vec[0];
    for (int i = 1; i < size; i++)
        if (res < vec[i]) res = vec[i];
    return res;
}

int countNegRadix(double num) {
    std::string str_num = std::to_string(num);
    if (str_num.find('.')) {
        while (str_num[str_num.size() - 1] == '0')  str_num.erase(str_num.size() - 1);
        return static_cast<int>(str_num.size()) - static_cast<int>(str_num.find('.')) - 1;
    } else {
        return 0;
    }
}

int getRemainder(double num) {
    int rad = countNegRadix(num);
    if (rad == 0) return 0;
    return static_cast<int>(round((((num - static_cast<int>(num)) * pow(10, rad)))));
}

std::vector<double> radixSortNeg(std::vector<double> vec, int radix) {
    std::vector<int> counts(10, 0);
    std::vector<int> counts1(10, 0);
    std::vector<double> res(vec.size());
    int size = static_cast<int>(vec.size());
    for (int i = 0; i < size; i++) {
        int r = getRemainder(vec[i]);
        counts[static_cast<int>((((r) % static_cast<unsigned int>(pow(10, (radix)))) / pow(10, (radix - 1))))]++;
    }
    for (int i = 1; i < 10; i++)
         for (int j = i - 1; j >= 0; j--)
             counts1[i] += counts[j];
    for (int i = 0; i < size; i++) {
        int r = getRemainder(vec[i]);
        int index = static_cast<int>((((r) % static_cast<unsigned int>(pow(10, (radix)))) / pow(10, (radix - 1))));
        res[counts1[index]] = vec[i];
        counts1[index]++;
    }
    return res;
}

std::vector<double> radixSortPos(std::vector<double> vec, int radix) {
    std::vector<int> counts(10, 0);
    std::vector<int> counts1(10, 0);
    std::vector<double> res(vec.size());
    int size = static_cast<int>(vec.size());
    for (int i = 0; i < size; i++) {
        int step = (static_cast<int>(vec[i]) % static_cast<unsigned int>(pow(10, (radix))));
        counts[static_cast<int>((step / pow(10, (radix - 1))))]++;
    }
    for (int i = 1; i < 10; i++)
        for (int j = i - 1; j >= 0; j--)
            counts1[i] += counts[j];
    for (int i = 0; i < size; i++) {
        int step = (static_cast<int>(vec[i]) % static_cast<unsigned int>(pow(10, (radix))));
        int index = static_cast<int>((step / pow(10, (radix - 1))));
        res[counts1[index]] = vec[i];
        counts1[index]++;
    }
    return res;
}

std::vector<double> radixSort(std::vector<double> vec) {
    int size = static_cast<int>(vec.size());
    int maxRadNeg = 0;
    int maxRadPos = countRadix(maxVal(vec));
    std::vector<double> res;
    for (int i = 0; i < size; i++)
        if (countNegRadix(vec[i]) > maxRadNeg) maxRadNeg = countNegRadix(vec[i]);
    if (maxRadNeg != 0) {
        res = radixSortNeg(vec, 1);
        if (maxRadNeg > 1) {
            for (int j = 2; j <= maxRadNeg; j++)
                res = radixSortNeg(res, j);
        }
        for (int radix = 1; radix <= maxRadPos; radix++)
            res = radixSortPos(res, radix);
    } else {
        res = radixSortPos(vec, 1);
        for (int radix = 2; radix <= maxRadPos; radix++)
            res = radixSortPos(res, radix);
    }
    return res;
}

void makeNetwork(int size) {
    cout << "makeNetWork in begin = " << endl;
    std::vector<int> vec(size);
    for (int i = 0; i < size; i++) {
        vec[i] = i;
    }
    cout << "makeNetWork before net = "<< endl;
    net(vec);
    cout << "makeNetWork after net = " << endl;
}

void net(std::vector<int> vec) {
    int size = static_cast<int>(vec.size());
    if (size <= 1) {
        return;
    }
    std::vector<int> left(vec.begin(), vec.begin() + size / 2);
    std::vector<int> right(vec.begin() + size / 2, vec.end());
    net(left);
    net(right);
    cout << "net before OddEvenMerge = " << endl;
    cout << "leftSize = "<< left.size() << endl;
    for (int i = 0; i < left.size(); i++) {
        cout << "left[" << i << "]= " << left[i] << endl;
    }
    cout << "rightSize = " << right.size() << endl;
    for (int i = 0; i < right.size(); i++) {
        cout << "right[" << i << "]= " << right[i] << endl;
    }
    oddEvenMerge(left, right);
    cout << "net afer OddEvenMerge = " << endl;
}

void oddEvenMerge(std::vector<int> left, std::vector<int> right) {
    int size = static_cast<int>(left.size()) + static_cast<int>(right.size());
    cout << "oddEvenMerge - Size = " << size << endl;
    if (size <= 1) {
        return;
    }
    if (size == 2) {
        cout << "oddEvenMerge in bloc size ==2 = " << endl;
        comps.push_back(std::pair<int, int>(left[0], right[0]));
        cout << "comps[0].first =  " << comps[0].first << endl;
        cout << "comps[0].first =  " << comps[0].second << endl;
        return;
    }
    std::vector<int> left_odd;
    std::vector<int> left_even;
    std::vector<int> right_odd;
    std::vector<int> right_even;
    int size_left = static_cast<int>(left.size());
    for (int i = 0; i < size_left; i++) {
        if (i % 2 == 0) {
            left_even.push_back(left[i]);;
        } else {
            left_odd.push_back(left[i]);;
        }
    }
    int size_right = static_cast<int>(right.size());
    for (int i = 0; i < size_right; i++) {
        if (i % 2 == 0) {
            right_even.push_back(right[i]);;
        } else {
            right_odd.push_back(right[i]);;
        }
    }

    oddEvenMerge(left_odd, right_odd);
    oddEvenMerge(left_even, right_even);

    std::vector<int> res(size_left+ size_right);
    for (int i = 0; i < size_left; i++) {
        res[i] = left[i];
    }
    for (int i = 0; i < size_right; i++) {
        res[i] = right[i];
    }

    for (int i = 1; i + 1 < size; i += 2) {
        comps.push_back(std::pair<int, int>(left[i], right[i]));
    }
}


std::vector<double> radix_sort_batcher_omp(std::vector<double> vec, int num_threads) {
    cout << "num_threads = " << num_threads << endl;
    omp_set_num_threads(num_threads);
    cout << "num_threads2 = " << num_threads << endl;
    int size = static_cast<int>(vec.size());
    makeNetwork(num_threads);
    cout << "num_threads3 = " << num_threads << endl;
    int addition = 0;
    while (size % num_threads != 0) {
            addition++;
            vec.push_back(0.0);
    }
    cout << "addition = " << addition << endl;
    int newSize = static_cast<int>(vec.size());
    cout << "newSize = " << newSize << endl;
    int localSize = newSize / num_threads;
    cout << "localSize = " << localSize << endl;
    std::vector<double> localVec(localSize);
    std::vector<double> localVec1(localSize);
    std::vector<double> localVec2(localSize);
#pragma omp parallel private(localVec, localVec1, localVec2) shared(vec, localSize)
    {
        int tid = omp_get_thread_num();
        cout << "tid = " << tid << endl;
        localVec.assign(vec.begin() + localSize * tid, vec.begin() + localSize * (tid + 1));
        localVec = radixSort(localVec);
#pragma omp barrier
        if (tid == 0) {
            for (int i = 0; i < localSize; i++) {
                cout << "local[" << i << "]= " << localVec[i] << endl;
            }
        }
#pragma omp barrier
        if (tid == 1) {
            for (int i = 0; i < localSize; i++) {
                cout << "local[" << i << "]= " << localVec[i] << endl;
            }
        }
#pragma omp barrier
        if (tid == 2) {
            for (int i = 0; i < localSize; i++) {
                cout << "local[" << i << "]= " << localVec[i] << endl;
            }
        }
        int begin = static_cast<int>(localSize * tid);
        if (tid == 0) {
            cout << "tid = 0, begin = " << begin << endl;
        }
        if (tid == 1) {
            cout << "tid = 1, begin = " << begin << endl;
        }
        if (tid == 2) {
            cout << "tid = 2, begin = " << begin << endl;
        }
        for (int i = 0; i < localSize; i++) {
            vec[begin + i] = localVec[i];
        }
        int countPair = static_cast<int>(comps.size());
        cout << "countPair" << countPair << endl;
        std::vector<double> localVec1(localSize);
        if (tid == 0) {

            for (int i = 0; i < countPair; i++) {
                cout << "comps.first = " << comps[i].first << "comps.second = " << comps[i].second << endl;
            }
        }
        std::vector<double> localVec2(localSize);
        for (int i = 0; i < countPair; i++) {
#pragma omp barrier
#pragma omp critical
                if (tid == comps[i].first) {
                    double tmp = vec[localSize * tid];
                    int local_point = localSize * tid;
                    int nearby_point = localSize * comps[i].second;
                    int step = 0;
                    for (int j = localSize * tid; j < localSize * tid + localSize; j++, step++) {
                        for (int h = localSize * comps[i].second + step; h < localSize * comps[i].second + localSize; h++) {
                            if (vec[j] < vec[h]) {
                                tmp = vec[j];
                            }
                            else {
                                vec[j] = vec[nearby_point];
                                nearby_point++;
                            }
                        }
                     }       
                    cout << "tid in for after inner for" << tid << endl;
                }
                else if (tid == comps[i].second) {
                    int start = localSize - 1;
                    int local_point = localSize * tid + localSize;
                    int nearby_point = localSize * comps[i].first + localSize;
                    for (int j = start; j >= 0; j--) {
                        int local = vec[local_point];
                        int nearby = vec[nearby_point];
                        if (local > nearby) {
                            localVec2[j] = local;
                            local_point--;
                        }
                        else {
                            localVec2[j] = nearby;
                            nearby_point--;
                        }
                    }
                    cout << "tid in for after inner for" << tid << endl;
                }
            
#pragma omp barrier
#pragma omp critical
                if (tid == comps[i].first) {
                    cout << "before vec assign1"<< endl;

                    cout << "size localVec1 = " << localVec1.size() << endl;
                    int j = 0;
                    for (int i =  localSize * tid; i < localSize * tid + localVec1.size(); i++)
                    {
                        cout << "vec["<<"] = " << localVec1[j] << endl;
                        vec[i] = localVec1[j];
                        j++;
                    }
                    //vec.assign((vec.begin() + localSize * tid), localVec1.begin());
                    cout << "after vec assign1"<< endl;
                }
#pragma omp barrier
#pragma omp critical
                if (tid == comps[i].second) {
                    cout << "before vec assign2" <<endl;
                    cout << "size localVec2 = " << localVec2.size() << endl;
                    int j = 0;
                    for (int i =  localSize * tid; i < localSize * tid + localVec2.size(); i++)
                    {
                        cout << "vec[" << "] = " << localVec2[j] << endl;
                        vec[i] = localVec2[j];
                        j++;
                    }
                    cout << "after vec assign2" << endl;
                }
#pragma omp barrier
        }
    }
    cout << "the end" << endl;
    return vec;
}

