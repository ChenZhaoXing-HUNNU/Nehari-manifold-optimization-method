function need_str = Get_deci(a,n)
%����һ��С�����õ�����С�����nλ�������������ַ����� 
num = floor(a);
str = num2str(num);
len = length(str);
str1 = num2str(a);
len_a =length(str1);
if len_a-len>=n+1
    need_str = str1(1:(len+n+1));
else
    need_str = num2str(a);
end