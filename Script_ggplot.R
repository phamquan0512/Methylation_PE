#thiet lap thu muc lam viec ( working directory)/ xem thu muc lam viec hien tai
setwd("D:/QUAN/CAO HOC/LUAN VAN/METHYLATION/nm")
getwd()
# tải packages
install.packages("ggplot2")  #lưu ý: dau ngoac kep
# truy cap packages
library(readxl)  # lưu ý: khong dau ngoac kep
library(datasets)
library(ggplot2)
library(gcookbook)
library(dplyr)
library(pastecs) # packages cho phan tích phân phối chuẩn???
remove.packages() # remove package

head(mtcars)
mpg
diamonds
iris

ggplot(iris,aes(x=Sepal.Length,y=Petal.Length))+
  labs(title = "Relationship between Sepal and Petal Length")

ggplot(iris,aes(x=Sepal.Length,y=Petal.Length))+
  labs(title = "Relationship between Sepal and Petal Length")+
  geom_point()

ggplot(iris,aes(x=Sepal.Length,fill=factor(Species)))+labs(title = "Distribution of Sepal Length")+geom_histogram(binwidth = 1)+facet_wrap(~Species)
ggplot(iris,aes(x=Species))+labs(title = "Count of each species")+geom_bar()
ggplot(iris,aes(x=Species,y=Sepal.Length,fill=Species))+labs(title = "Sepal length of each species")+geom_boxplot()



# để lọc những dữ liệu mình cần, có thể xài filter
#vd: trong data diamonds, lọc những dòng có cut là ideal
ideal=filter(diamonds,cut=="Ideal")
#vd: trong data diamonds, lọc những dòng có cut là Ideal và Fair, color là E và I, price > 15000
ideal=filter(diamonds,cut==c("Ideal","Fair"),color==c("E","I"),price >15000)


baseplot<-ggplot(mtcars, aes(x=hp,y=cyl))+labs(title = "HP~Cylon", x="HP", y="Cylon")
print(baseplot)
finalplot<-baseplot + geom_function()
#barplot
baseplot<-ggplot(mtcars, aes(x=as.factor(cyl)))+labs(title = "Number of car using different type of cylinder", x="Cylinder type", y="Number of car")+geom_bar()
print(baseplot)
#piechart
count(data,x) #x là biến muốn tính %
mutate(data,prop = n / sum(n) * 100) # tính % theo từng giá trị của biến và thêm vào thành một cột trong data
cyl_data <- mtcars %>%
  +     count(cyl) %>%
  +     mutate(prop = n / sum(n) * 100) 

ggplot(cyl_data, aes(x = "", y = prop, fill = factor(cyl))) +
  geom_bar(stat = "identity", width = 1) +
  coord_polar(theta = "y") +
  labs(fill = "Cylinders", title = "Proportion of Cylinders in mtcars Dataset") +
  theme_void()
#homework: dùng diamond data, hãy cho biết % của mỗi loại cut
cut_data<-diamonds %>%
count(cut) %>%
mutate(prop = n/sum(n)*100)
cut_data
ggplot(cut_data,aes(x="",y=prop,fill=cut))+geom_bar(stat = "identity",width=1)+coord_polar(theta = "y")+labs(title="Proportion of diamond cut")+theme_void()
#homework: hãy cho biết % của mỗi loại clarity
clarity_data<-diamonds %>%
  count(clarity) %>%
  mutate(prop=n/sum(n)*100)
clarity_data
ggplot(clarity_data,aes(x="",y=prop,fill=clarity))+geom_bar(stat = "identity",width = 1)+coord_polar(theta = "y")+labs(title="Proportion of diamond's clarity")+theme_void()

## histogram
baseplot<-ggplot(mtcars,aes(x=biến liên tục))+geom_histogram(binwidth = )+labs(title="...")
print(baseplot)
b<-diamonds
head(b)
ggplot(b,aes(x=price))+geom_histogram(binwidth = 5)+labs(title="Distribution of diamond price")
ggplot(b,aes(x=depth))+geom_histogram(binwidth = 1)+labs(title="Distribution of diamond depth")
ggplot(b,aes(x=price,fill=factor(cut)))+geom_histogram(binwidth = 5)+labs(title="Distribution of diamond price")
#chức năng facet_wrap(~cut)/ facet_grid(~cut): tách các giá trị price theo cut thành từng cut riêng
ggplot(b,aes(x=price,fill=factor(cut)))+geom_histogram(binwidth = 50)+labs(title="Distribution of diamond price")+facet_wrap(~cut)

## density plot
# đánh giá s
baseplot<-ggplot(data,aes(x=biến liên tục))+geom_density()+labs(title="...")
print(baseplot)
# visualize distribution of diamond's price based on clarity, color and cut
ggplot(b,aes(x=price,fill=factor(clarity)))+geom_density()+labs(title="Distribution of diamond's price based on clarity")
ggplot(b,aes(x=price,fill=factor(color)))+geom_density()+labs(title="Distribution of diamond's price based on color")
ggplot(b,aes(x=price,fill=factor(cut)))+geom_density()+labs(title="Distribution of diamond's price based on cut")

# Create the heatmap
ggplot(mtcars_melted, aes(x = Variable, y = Car, fill = Value))+geom_tile(color = "white")+scale_fill_gradient(low = "green", high = "red")+theme_minimal()+theme(axis.text.x = element_text(angle = 45, hjust = 1))+labs(title = "Heatmap of mtcars Dataset", x = "Variable", y = "Car")
savehistory("D:/QUAN/CAO HOC/Microbiome_course/GGPLOT2/ggplot.Rhistory")
