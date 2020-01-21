ggplot()+
  geom_point(aes(x=rep(1,10), y=1:10)) + geom_line(aes(x=rep(1,10), y=1:10))+
  geom_point(aes(x=rep(2,10), y=1:10)) + geom_line(aes(x=rep(2,10), y=1:10))+
  geom_point(aes(x=rep(3,10), y=1:10)) + geom_line(aes(x=rep(3,10), y=1:10))+
  geom_point(aes(x=rep(4,10), y=1:10)) + geom_line(aes(x=rep(4,10), y=1:10))+
  geom_point(aes(x=rep(5,10), y=1:10)) + geom_line(aes(x=rep(5,10), y=1:10))+
  geom_point(aes(x=rep(6,10), y=1:10)) + geom_line(aes(x=rep(6,10), y=1:10))+
  geom_point(aes(x=rep(7,10), y=1:10)) + geom_line(aes(x=rep(7,10), y=1:10))+
  geom_point(aes(x=rep(8,10), y=1:10)) + geom_line(aes(x=rep(8,10), y=1:10))+
  geom_point(aes(x=rep(9,10), y=1:10)) + geom_line(aes(x=rep(9,10), y=1:10))+
  geom_point(aes(x=rep(10,10), y=1:10)) + geom_line(aes(x=rep(10,10), y=1:10))+
  geom_point(aes(x=1:10, y=rep(1,10))) + geom_line(aes(x=1:10, y=rep(1,10)))+
  geom_point(aes(x=1:10, y=rep(2,10))) + geom_line(aes(x=1:10, y=rep(2,10)))+
  geom_point(aes(x=1:10, y=rep(3,10))) + geom_line(aes(x=1:10, y=rep(3,10)))+
  geom_point(aes(x=1:10, y=rep(4,10))) + geom_line(aes(x=1:10, y=rep(4,10)))+
  geom_point(aes(x=1:10, y=rep(5,10))) + geom_line(aes(x=1:10, y=rep(5,10)))+
  geom_point(aes(x=1:10, y=rep(6,10))) + geom_line(aes(x=1:10, y=rep(6,10)))+
  geom_point(aes(x=1:10, y=rep(7,10))) + geom_line(aes(x=1:10, y=rep(7,10)))+
  geom_point(aes(x=1:10, y=rep(8,10))) + geom_line(aes(x=1:10, y=rep(8,10)))+
  geom_point(aes(x=1:10, y=rep(9,10))) + geom_line(aes(x=1:10, y=rep(9,10)))+
  geom_point(aes(x=1:10, y=rep(10,10))) + geom_line(aes(x=1:10, y=rep(10,10)))+
  xlab("x") + ylab("y")