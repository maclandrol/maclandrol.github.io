---
layout: page
title: Blog posts
---

<div>
<ul class="post-list">
{% for post in site.posts limit:10 %} 
	
	<li>{{ post.date | date_to_string }} &raquo; <a href="{{ post.url }}">{{ post.title }}</a>
	</li>
{% endfor %}
	</ul>

</div>

<div class="infinite-spinner"></div>