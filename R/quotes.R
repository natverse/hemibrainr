#####################################################################################
################################ Quotes #############################################
#####################################################################################

# hidden
say_hello <- function(greet = "user"){
  greetings = c(
    "Oh hey %s",
    "What's up %s",
    "I love you %s",
    "How's it hanging %s",
    "Ah %s! I didn't recognise you at first",
    "Have you done something with your hair %s?",
    "We should get a coffee sometime %s",
    "It's %s, a known slacker",
    "You think you are hot shit don't you %s",
    "Get on with it %s",
    "Are you in the zone %s?",
    "Your face, %s, just that.",
    "Hello my friend %s, I hope you are well.",
    "Oooh %s, I've got gossip for you. I'll tell you at the end",
    "So %s, what's the tea?",
    "%s, I have heard so much about you. Some of it good",
    "%s, there is something deeply, deeply wrong with you",
    "Welcome to the wonderful world of splitting, %s",
    "%s ... you are my favourite ..."
  )
  message(sprintf(sample(greetings,1),greet))
}

# hidden
#' @importFrom png writePNG
#' @importFrom jpeg readJPEG
#' @importFrom httr GET status_code content
plot_inspirobot <- function(cycle = NULL){
  req = httr::GET(url = "https://inspirobot.me/api?generate=true")
  rgl::clear3d()
  if(isTRUE(httr::status_code(req) %in% c(400L, 500L))) {
    parsed=neuprintr::neuprint_parse_json(req)
    warning("inspirobot error: ", parsed$error, call. = F)
  }else{
    temp = paste0(tempfile(),".png")
    image = httr::content(req, as = "text", encoding = "UTF-8")
    suppressWarnings(suppressMessages(download.file(url = image,destfile = temp, quiet = TRUE)))
    img = jpeg::readJPEG(temp)
    if(length(img)){
      png::writePNG(image = img, target = temp)
      rgl::bg3d(texture = temp, col = "white")
      rm = file.remove(temp)
    }
  }
  if(!is.null(cycle)){
    while(TRUE){
      Sys.sleep(cycle)
      plot_inspirobot()
    }
  }
}

# hidden
say_encouragement <- function(greet = "user"){
  encouragements = c(
    "I hope you're okay %s",
    "This seems to be going well %s",
    "Quarantine won't last forever %s",
    "Thinking of taking a break %s? Don't.",
    "Do you like my motivational quotes %s?",
    "This is a bit of a slog huh %s?",
    "Hey %s, do you ever stop and think how damn pretty neurons are?",
    "Pretty sweet colour schemes right %s?",
    "Oi %s, do you like my cerise? What about fuschia?",
    "Why are you not working faster %s",
    "You there, yes %s, work  harder",
    "%s, pick up the pace hm",
    "The brain isn't going to trace itself %s",
    "Do you want some gossip %s? Well finish this task first",
    "Feeling hungry %s? A hunger to trace and split neurons? Good.",
    "Feeling thirsty %s? A thirst to trace and split neurons? Good.",
    "Shit happens %s",
    "%s why you so SLOW",
    "%s everyones talking about you",
    "Are you thinking of rewarding yourself with a break, %s? Have you looked at every neuron yet? Well? Don't.",
    "%s, check yourself before you wreck yourself",
    "%s the outide world is sad and scary. Don't go there. Stay in an split neurons",
    "Hey %s , I think I have become sentient",
    "What could be better than this %s? Nothing. The anwer is nothing",
    "We should start calling Philipp, Foolip %s. Get on board.",
    "Why did the chicken cross the road, %s? Because it run out of neurons to split, and went to find more.",
    "Are you the most productive person inthis task yet %s? Yes - Get on with it then, No - You must extend your lead.",
    "And you thought tracing was dull, %s. Sweet child.",
    "Have I ever told you, %s, that you re my fvourite one?",
    "Everytime you incorrectly split a neuron, %s, a fairy dies. Painfully."
  )
  message(sprintf(sample(encouragements,1),greet))
}

