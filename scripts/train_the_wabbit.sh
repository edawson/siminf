cat $1 | shuf |  vw --passes=20 --ect 5 --cache_file cache.f --binary --interactions vvvv -f trained.model
