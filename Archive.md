---
layout: page
title: Blog posts
---

{% if paginator.page %}
   {% assign offset = paginator.page | minus:1 | times:paginator.per_page %} 
   {% assign currentPage = paginator.page %}
{% else %} 
  {% assign offset = 0 %} 
  {% assign currentPage = 1 %}
{% endif %}

<div id="blogContainer" data-page="{{ currentPage }}" data-totalPages="{{ paginator.total_pages }}">

{% for post in site.posts limit:site.paginate offset:offset %} 
  <ul>
  <li>{{ post.date | date_to_string }} &raquo; <a href="{{ post.url }}">{{ post.title }}</a>
  </li>
  </ul>
{% endfor %}


</div>

{% assign postCount = site.posts | size %} 
{% assign postsCovered = site.paginate | plus:offset %}

<div class="pagination">
   {% if postsCovered < postCount %}
     <a class="pagination-item newer">Load more</a>
  {% endif %}
</div>

<script>
  {% include pagination.js %}
</script>
