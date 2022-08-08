# On terminal (after installing docker):
# sudo systemctl start docker
# sudo systemctl enable docker
# sudo docker run -d -p 4445:4444 selenium/standalone-chrome

library(RSelenium)
library(tidyverse)
library(rvest)
library(here)
remDr <- remoteDriver(port=4445L, browserName="chrome")
remDr$open()

# Connect to website
remDr$navigate("https://www.rekihaku.ac.jp/up-cgi/login.pl?p=param/jomo/db_param")
# Choose 500 from the dropdown list
webElem <- remDr$findElement(using = 'xpath', value = "/html/body/form/table[2]/tbody/tr/td[2]/select/option[3]")
webElem$clickElement()

# Blank field search to retrieve all sites
address_element <- remDr$findElement(using = 'css', value = 'tr:nth-child(8) input')
address_element$sendKeysToElement(list(" "))
button_element <- remDr$findElement(using = 'css', value = "center input:nth-child(1)")
button_element$clickElement()

# Recover total number of hits
session=read_html(remDr$getPageSource()[[1]])
nhits=html_text(html_nodes(session,'*')[4])
nhits=as.numeric(unlist(strsplit(strsplit(nhits,'検索結果：')[[1]][2],'件データ')[[1]][1]))
nclicks=ceiling(nhits/500)
# Recover table from the first page
table <- session %>% html_nodes("table")
tmp.out.df=html_table(table[3])[[1]] #covers 1-50

# Recover tables from the other pages
  for (i in 2:nclicks)
  {
    Sys.sleep(3)
    print(i)
    if (i==2) {  button_element <- remDr$findElement(using = 'css', value = "input:nth-child(9)")}
    if (i>2) {  button_element <- remDr$findElement(using = 'css', value = "input:nth-child(10)")}
    button_element$clickElement()
    session=read_html(remDr$getPageSource()[[1]])
    table <- session %>% html_nodes("table")
    tmp.out.df=rbind.data.frame(tmp.out.df,html_table(table[3])[[1]])
  }

write.csv(tmp.out.df,file=here('data','site_raw.csv'),row.names=FALSE)
