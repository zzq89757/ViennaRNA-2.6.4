#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// 定义结构体
typedef struct {
    int age;
    double height;
    char name[50];
} Person;

void personToString(duplexT mfe, char *result, size_t resultSize) {
    snprintf(result, resultSize, "%d|%f|%s", person.age, person.height, person.name);
}

int main() {
    // 创建一个结构体实例
    Person p;
    p.age = 30;
    p.height = 5.9;
    strcpy(p.name, "John Doe");
    
    // 存储结果字符串
    char result[100];
    
    // 调用函数，将结构体属性转换为分隔的字符串
    personToString(p, result, sizeof(result));
    
    // 打印结果
    printf("Resulting string: %s\n", result);
    
    return 0;
}